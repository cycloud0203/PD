// ObjectiveFunction.cpp
#include "ObjectiveFunction.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>

namespace {

constexpr double kSmoothingGamma = 0.1;
constexpr double kMinExpSum = 1e-10;
constexpr double kMinLambdaMagnitude = 1e-20;
constexpr double kDefaultDensityTargetFactor = 0.8;
constexpr double kMinContribution = 1e-10;

int resolveGridCount(int gridCount, unsigned numModules) {
    if (gridCount > 0) {
        return gridCount;
    }
    return static_cast<int>(std::sqrt(static_cast<double>(numModules)));
}

double bellOverlap(double distance, double elementSize, double unitSize) {
    const double nearLimit = elementSize * 0.5 + unitSize;
    const double farLimit = elementSize * 0.5 + 2.0 * unitSize;
    const double denomNear = elementSize + 2.0 * unitSize;
    const double denomFar = elementSize + 4.0 * unitSize;
    const double quadCoef = 4.0 / (denomNear * denomFar);
    const double linCoef = 2.0 / (unitSize * denomFar);

    if (distance <= nearLimit) {
        return 1.0 - quadCoef * distance * distance;
    }
    if (distance <= farLimit) {
        return linCoef * (distance - farLimit);
    }
    return 0.0;
}

}  // namespace

Wirelength::Wirelength(Placement& placement)
    : BaseFunction(placement.numModules()),
      placement_(&placement),
      gamma_(kSmoothingGamma) {}

double Wirelength::calculate() {
    const int netCount = static_cast<int>(placement_->numNets());
    double totalLength = 0.0;
    const double gamma = gamma_;

    auto smoothMax = [&](const std::vector<double>& values) {
        const double maxValue = *std::max_element(values.begin(), values.end());
        double weightSum = 0.0;
        double weightedSum = 0.0;
        for (double value : values) {
            const double weight = std::exp((value - maxValue) / gamma);
            weightSum += weight;
            weightedSum += value * weight;
        }
        return weightedSum / weightSum;
    };

    auto smoothMin = [&](const std::vector<double>& values) {
        const double minValue = *std::min_element(values.begin(), values.end());
        double weightSum = 0.0;
        double weightedSum = 0.0;
        for (double value : values) {
            const double weight = std::exp((minValue - value) / gamma);
            weightSum += weight;
            weightedSum += value * weight;
        }
        return weightedSum / weightSum;
    };

#pragma omp parallel for reduction(+ : totalLength) schedule(static)
    for (int netId = 0; netId < netCount; ++netId) {
        Net& net = placement_->net(netId);
        const int pinCount = static_cast<int>(net.numPins());
        if (pinCount < 2) {
            continue;
        }

        std::vector<double> xs;
        std::vector<double> ys;
        xs.reserve(pinCount);
        ys.reserve(pinCount);
        for (int pinIndex = 0; pinIndex < pinCount; ++pinIndex) {
            Pin& pin = net.pin(pinIndex);
            xs.push_back(pin.x());
            ys.push_back(pin.y());
        }

        totalLength += (smoothMax(xs) - smoothMin(xs)) + (smoothMax(ys) - smoothMin(ys));
    }

    value_ = totalLength;
    return value_;
}

void Wirelength::Backward() {
    const int moduleCount = static_cast<int>(placement_->numModules());
    const int netCount = static_cast<int>(placement_->numNets());

    std::fill(grad_.begin(), grad_.end(), Point2<double>{0.0, 0.0});
    value_ = 0.0;

#pragma omp parallel
    {
        std::vector<Point2<double>> localForces(moduleCount, Point2<double>{0.0, 0.0});
        double localWirelength = 0.0;

#pragma omp for schedule(static)
        for (int netIndex = 0; netIndex < netCount; ++netIndex) {
            Net& net = placement_->net(netIndex);
            const int pinCount = static_cast<int>(net.numPins());
            if (pinCount <= 1) {
                continue;
            }

            double right = net.pin(0).x();
            double left = right;
            double top = net.pin(0).y();
            double bottom = top;

            for (int pinIndex = 1; pinIndex < pinCount; ++pinIndex) {
                Pin& pin = net.pin(pinIndex);
                right = std::max(right, pin.x());
                left = std::min(left, pin.x());
                top = std::max(top, pin.y());
                bottom = std::min(bottom, pin.y());
            }

            const double invGamma = 1.0 / gamma_;
            double weightSums[4] = {0.0, 0.0, 0.0, 0.0};
            double expSums[4] = {0.0, 0.0, 0.0, 0.0};

            for (int pinIndex = 0; pinIndex < pinCount; ++pinIndex) {
                Pin& pin = net.pin(pinIndex);
                const double x = pin.x();
                const double y = pin.y();

                const double expValues[4] = {
                    std::exp((x - right) * invGamma),
                    std::exp((left - x) * invGamma),
                    std::exp((y - top) * invGamma),
                    std::exp((bottom - y) * invGamma)};
                const double posValues[4] = {x, -x, y, -y};

                for (int axis = 0; axis < 4; ++axis) {
                    weightSums[axis] += posValues[axis] * expValues[axis];
                    expSums[axis] += expValues[axis];
                }
            }

            double netWirelength = 0.0;
            for (int axis = 0; axis < 4; ++axis) {
                if (expSums[axis] > kMinExpSum) {
                    netWirelength += weightSums[axis] / expSums[axis];
                }
            }
            localWirelength += netWirelength;

            for (int pinIndex = 0; pinIndex < pinCount; ++pinIndex) {
                Pin& pin = net.pin(pinIndex);
                const int moduleId = static_cast<int>(pin.moduleId());
                if (placement_->module(moduleId).isFixed()) {
                    continue;
                }

                const double x = pin.x();
                const double y = pin.y();
                const double expValues[4] = {
                    std::exp((x - right) * invGamma),
                    std::exp((left - x) * invGamma),
                    std::exp((y - top) * invGamma),
                    std::exp((bottom - y) * invGamma)};
                const double posValues[4] = {x, -x, y, -y};
                double commonTerms[4];
                for (int axis = 0; axis < 4; ++axis) {
                    commonTerms[axis] = expValues[axis] * invGamma;
                }

                if (expSums[0] > kMinExpSum) {
                    const double numerator =
                        expSums[0] * (expValues[0] + posValues[0] * commonTerms[0]);
                    const double denominator = weightSums[0] * commonTerms[0];
                    localForces[moduleId].x +=
                        (numerator - denominator) / (expSums[0] * expSums[0]);
                }
                if (expSums[1] > kMinExpSum) {
                    const double numerator =
                        expSums[1] * (expValues[1] + posValues[1] * commonTerms[1]);
                    const double denominator = weightSums[1] * commonTerms[1];
                    localForces[moduleId].x -=
                        (numerator - denominator) / (expSums[1] * expSums[1]);
                }
                if (expSums[2] > kMinExpSum) {
                    const double numerator =
                        expSums[2] * (expValues[2] + posValues[2] * commonTerms[2]);
                    const double denominator = weightSums[2] * commonTerms[2];
                    localForces[moduleId].y +=
                        (numerator - denominator) / (expSums[2] * expSums[2]);
                }
                if (expSums[3] > kMinExpSum) {
                    const double numerator =
                        expSums[3] * (expValues[3] + posValues[3] * commonTerms[3]);
                    const double denominator = weightSums[3] * commonTerms[3];
                    localForces[moduleId].y -=
                        (numerator - denominator) / (expSums[3] * expSums[3]);
                }
            }
        }

#pragma omp critical
        {
            for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
                grad_[moduleId].x += localForces[moduleId].x;
                grad_[moduleId].y += localForces[moduleId].y;
            }
            value_ += localWirelength;
        }
    }
}

Density::Density(Placement& placement, int* gridCountPtr)
    : BaseFunction(placement.numModules()),
      placement_(&placement),
      gridCountPtr_(gridCountPtr),
      targetFactor_(kDefaultDensityTargetFactor) {}

double Density::calculate() {
    const int gridSize =
        resolveGridCount(gridCountPtr_ ? *gridCountPtr_ : 0, placement_->numModules());

    const double minX = placement_->boundryLeft();
    const double maxX = placement_->boundryRight();
    const double minY = placement_->boundryBottom();
    const double maxY = placement_->boundryTop();

    const double regionWidth = maxX - minX;
    const double regionHeight = maxY - minY;
    const double cellWidth = regionWidth / gridSize;
    const double cellHeight = regionHeight / gridSize;
    const double capacityLimit = targetFactor_ * cellWidth * cellHeight;
    const double cellArea = cellWidth * cellHeight;

    const int cellCount = gridSize * gridSize;
    std::vector<double> gridDensity(cellCount, 0.0);

#pragma omp parallel
    {
        std::vector<double> localDensity(cellCount, 0.0);

#pragma omp for schedule(static)
        for (int moduleId = 0; moduleId < static_cast<int>(placement_->numModules()); ++moduleId) {
            Module& module = placement_->module(moduleId);
            if (module.isFixed()) {
                continue;
            }

            const double centerX = module.centerX();
            const double centerY = module.centerY();
            const double width = module.width();
            const double height = module.height();

            const int colStart = std::max(
                0, static_cast<int>((centerX - width * 0.5 - minX) / cellWidth) - 2);
            const int colEnd = std::min(
                gridSize - 1,
                static_cast<int>((centerX + width * 0.5 - minX) / cellWidth) + 2);
            const int rowStart = std::max(
                0, static_cast<int>((centerY - height * 0.5 - minY) / cellHeight) - 2);
            const int rowEnd = std::min(
                gridSize - 1,
                static_cast<int>((centerY + height * 0.5 - minY) / cellHeight) + 2);

            for (int row = rowStart; row <= rowEnd; ++row) {
                const double cellCenterY = minY + (row + 0.5) * cellHeight;
                const double verticalDistance = std::abs(centerY - cellCenterY);
                const double overlapY =
                    bellOverlap(verticalDistance, height, cellHeight);

                for (int col = colStart; col <= colEnd; ++col) {
                    const double cellCenterX = minX + (col + 0.5) * cellWidth;
                    const double horizontalDistance = std::abs(centerX - cellCenterX);
                    const double overlapX =
                        bellOverlap(horizontalDistance, width, cellWidth);
                    const double contribution = cellArea * overlapX * overlapY;
                    if (contribution > kMinContribution) {
                        localDensity[row * gridSize + col] += contribution;
                    }
                }
            }
        }

#pragma omp critical
        {
            for (int index = 0; index < cellCount; ++index) {
                gridDensity[index] += localDensity[index];
            }
        }
    }

    double totalPenalty = 0.0;
#pragma omp parallel for reduction(+ : totalPenalty) schedule(static)
    for (int index = 0; index < cellCount; ++index) {
        const double excess = gridDensity[index] - capacityLimit;
        if (excess > 0.0) {
            totalPenalty += excess * excess;
        }
    }

    value_ = totalPenalty;
    return value_;
}

void Density::Backward() {
    const int moduleCount = static_cast<int>(placement_->numModules());
    std::fill(grad_.begin(), grad_.end(), Point2<double>{0.0, 0.0});

    const int gridSize =
        resolveGridCount(gridCountPtr_ ? *gridCountPtr_ : 0, placement_->numModules());

    const double left = placement_->boundryLeft();
    const double right = placement_->boundryRight();
    const double bottom = placement_->boundryBottom();
    const double top = placement_->boundryTop();
    const double regionWidth = right - left;
    const double regionHeight = top - bottom;
    const double cellWidth = regionWidth / gridSize;
    const double cellHeight = regionHeight / gridSize;
    const double capacityLimit = targetFactor_ * cellWidth * cellHeight;
    const double cellArea = cellWidth * cellHeight;

    std::vector<std::vector<double>> gridDensity(
        gridSize, std::vector<double>(gridSize, 0.0));

    const int threadCount = omp_get_max_threads();
    std::vector<std::vector<std::vector<double>>> threadDensityContributions(threadCount);
    for (int threadId = 0; threadId < threadCount; ++threadId) {
        threadDensityContributions[threadId].assign(
            gridSize, std::vector<double>(gridSize, 0.0));
    }

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        std::vector<double> cellContribution(gridSize * gridSize, 0.0);

#pragma omp for schedule(static)
        for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
            if (placement_->module(moduleId).isFixed()) {
                continue;
            }

            const double centerX = placement_->module(moduleId).centerX();
            const double centerY = placement_->module(moduleId).centerY();
            const double width = placement_->module(moduleId).width();
            const double height = placement_->module(moduleId).height();
            const double halfWidth = width * 0.5;
            const double halfHeight = height * 0.5;

            const double invCellWidth = 1.0 / cellWidth;
            const double invCellHeight = 1.0 / cellHeight;

            const int xStart = std::max(
                0, static_cast<int>(std::floor((centerX - halfWidth - left) * invCellWidth)) - 2);
            const int xEnd = std::min(
                gridSize - 1,
                static_cast<int>(std::ceil((centerX + halfWidth - left) * invCellWidth)) + 2);
            const int yStart = std::max(
                0,
                static_cast<int>(std::floor((centerY - halfHeight - bottom) * invCellHeight)) - 2);
            const int yEnd = std::min(
                gridSize - 1,
                static_cast<int>(std::ceil((centerY + halfHeight - bottom) * invCellHeight)) + 2);

            const double widthNear = halfWidth + cellWidth;
            const double widthFar = halfWidth + 2.0 * cellWidth;
            const double heightNear = halfHeight + cellHeight;
            const double heightFar = halfHeight + 2.0 * cellHeight;

            const double widthDenomNear = width + 2.0 * cellWidth;
            const double widthDenomFar = width + 4.0 * cellWidth;
            const double densityAX = 4.0 / (widthDenomNear * widthDenomFar);
            const double densityBX = 2.0 / (cellWidth * widthDenomFar);

            const double heightDenomNear = height + 2.0 * cellHeight;
            const double heightDenomFar = height + 4.0 * cellHeight;
            const double densityAY = 4.0 / (heightDenomNear * heightDenomFar);
            const double densityBY = 2.0 / (cellHeight * heightDenomFar);

            constexpr int kBlockSize = 4;
            for (int blockX = xStart; blockX <= xEnd; blockX += kBlockSize) {
                const int blockXEnd = std::min(blockX + kBlockSize - 1, xEnd);
                for (int blockY = yStart; blockY <= yEnd; blockY += kBlockSize) {
                    const int blockYEnd = std::min(blockY + kBlockSize - 1, yEnd);
                    for (int gridX = blockX; gridX <= blockXEnd; ++gridX) {
                        const double gridCenterX = left + (gridX + 0.5) * cellWidth;
                        const double deltaX = std::abs(centerX - gridCenterX);
                        if (deltaX > widthFar) {
                            continue;
                        }

                        double overlapX = 0.0;
                        if (deltaX <= widthNear) {
                            overlapX = 1.0 - densityAX * deltaX * deltaX;
                        } else {
                            overlapX = densityBX * (deltaX - widthFar);
                        }
                        if (overlapX <= 0.0) {
                            continue;
                        }

                        for (int gridY = blockY; gridY <= blockYEnd; ++gridY) {
                            const double gridCenterY = bottom + (gridY + 0.5) * cellHeight;
                            const double deltaY = std::abs(centerY - gridCenterY);
                            if (deltaY > heightFar) {
                                continue;
                            }

                            double overlapY = 0.0;
                            if (deltaY <= heightNear) {
                                overlapY = 1.0 - densityAY * deltaY * deltaY;
                            } else {
                                overlapY = densityBY * (deltaY - heightFar);
                            }

                            if (overlapY > 0.0) {
                                cellContribution[gridX * gridSize + gridY] +=
                                    cellArea * overlapX * overlapY;
                            }
                        }
                    }
                }
            }
        }

        for (int gridX = 0; gridX < gridSize; ++gridX) {
            for (int gridY = 0; gridY < gridSize; ++gridY) {
                const double contribution = cellContribution[gridX * gridSize + gridY];
                if (contribution > 0.0) {
                    threadDensityContributions[threadId][gridX][gridY] += contribution;
                }
            }
        }
    }

    for (int threadId = 0; threadId < threadCount; ++threadId) {
        for (int gridX = 0; gridX < gridSize; ++gridX) {
            for (int gridY = 0; gridY < gridSize; ++gridY) {
                gridDensity[gridX][gridY] += threadDensityContributions[threadId][gridX][gridY];
            }
        }
    }

    std::vector<std::vector<Point2<double>>> threadGradients(threadCount);
    for (int threadId = 0; threadId < threadCount; ++threadId) {
        threadGradients[threadId].assign(moduleCount, Point2<double>{0.0, 0.0});
    }

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        std::vector<double> localGradX(moduleCount, 0.0);
        std::vector<double> localGradY(moduleCount, 0.0);

#pragma omp for schedule(static)
        for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
            if (placement_->module(moduleId).isFixed()) {
                continue;
            }

            const double centerX = placement_->module(moduleId).centerX();
            const double centerY = placement_->module(moduleId).centerY();
            const double width = placement_->module(moduleId).width();
            const double height = placement_->module(moduleId).height();
            const double halfWidth = width * 0.5;
            const double halfHeight = height * 0.5;

            const double invCellWidth = 1.0 / cellWidth;
            const double invCellHeight = 1.0 / cellHeight;

            const double widthDenomNear = width + 2.0 * cellWidth;
            const double widthDenomFar = width + 4.0 * cellWidth;
            const double quadCoefX = 4.0 / (widthDenomNear * widthDenomFar);
            const double linCoefX = 2.0 / (cellWidth * widthDenomFar);

            const double heightDenomNear = height + 2.0 * cellHeight;
            const double heightDenomFar = height + 4.0 * cellHeight;
            const double quadCoefY = 4.0 / (heightDenomNear * heightDenomFar);
            const double linCoefY = 2.0 / (cellHeight * heightDenomFar);

            const int xMin = std::max(
                0, static_cast<int>((centerX - halfWidth - left) * invCellWidth) - 2);
            const int xMax = std::min(
                gridSize - 1,
                static_cast<int>((centerX + halfWidth - left) * invCellWidth) + 2);
            const int yMin = std::max(
                0, static_cast<int>((centerY - halfHeight - bottom) * invCellHeight) - 2);
            const int yMax = std::min(
                gridSize - 1,
                static_cast<int>((centerY + halfHeight - bottom) * invCellHeight) + 2);

            const double radiusXNear = halfWidth + cellWidth;
            const double radiusXFar = halfWidth + 2.0 * cellWidth;
            const double radiusYNear = halfHeight + cellHeight;
            const double radiusYFar = halfHeight + 2.0 * cellHeight;
            const double cellFactor = 2.0 * cellArea;

            constexpr int kBlockSize = 2;
            double totalGradX = 0.0;
            double totalGradY = 0.0;

            for (int blockX = xMin; blockX <= xMax; blockX += kBlockSize) {
                for (int blockY = yMin; blockY <= yMax; blockY += kBlockSize) {
                    for (int xOffset = 0; xOffset < kBlockSize && blockX + xOffset <= xMax; ++xOffset) {
                        const int gridX = blockX + xOffset;
                        const double cellX = left + (gridX + 0.5) * cellWidth;
                        const double deltaX = std::abs(centerX - cellX);
                        if (deltaX > radiusXFar) {
                            continue;
                        }

                        double overlapX = 0.0;
                        double gradX = 0.0;
                        if (deltaX <= radiusXNear) {
                            overlapX = 1.0 - quadCoefX * deltaX * deltaX;
                            gradX = (centerX >= cellX) ? -2.0 * quadCoefX * deltaX
                                                         : 2.0 * quadCoefX * deltaX;
                        } else {
                            const double offsetX = deltaX - radiusXFar;
                            overlapX = linCoefX * offsetX;
                            gradX = (centerX >= cellX) ? 2.0 * linCoefX * offsetX
                                                       : -2.0 * linCoefX * offsetX;
                        }

                        for (int yOffset = 0; yOffset < kBlockSize && blockY + yOffset <= yMax;
                             ++yOffset) {
                            const int gridY = blockY + yOffset;
                            const double cellY = bottom + (gridY + 0.5) * cellHeight;
                            const double deltaY = std::abs(centerY - cellY);
                            if (deltaY > radiusYFar) {
                                continue;
                            }

                            double overlapY = 0.0;
                            double gradY = 0.0;
                            if (deltaY <= radiusYNear) {
                                overlapY = 1.0 - quadCoefY * deltaY * deltaY;
                                gradY = (centerY >= cellY) ? -2.0 * quadCoefY * deltaY
                                                           : 2.0 * quadCoefY * deltaY;
                            } else {
                                const double offsetY = deltaY - radiusYFar;
                                overlapY = linCoefY * offsetY;
                                gradY = (centerY >= cellY) ? 2.0 * linCoefY * offsetY
                                                           : -2.0 * linCoefY * offsetY;
                            }

                            const double overflow = gridDensity[gridX][gridY] - capacityLimit;
                            if (overflow > 0.0) {
                                const double force = cellFactor * overflow;
                                totalGradX += force * gradX * overlapY;
                                totalGradY += force * overlapX * gradY;
                            }
                        }
                    }
                }
            }

            localGradX[moduleId] = totalGradX;
            localGradY[moduleId] = totalGradY;
        }

        for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
            threadGradients[threadId][moduleId].x = localGradX[moduleId];
            threadGradients[threadId][moduleId].y = localGradY[moduleId];
        }
    }

    for (int threadId = 0; threadId < threadCount; ++threadId) {
        for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
            grad_[moduleId].x += threadGradients[threadId][moduleId].x;
            grad_[moduleId].y += threadGradients[threadId][moduleId].y;
        }
    }
}

ObjectiveFunction::ObjectiveFunction(Placement& placement)
    : BaseFunction(placement.numModules()),
      placement_(&placement),
      densityGridCount_(0),
      lambda_(1.0) {
    wirelength_ = std::make_unique<Wirelength>(placement);
    density_ = std::make_unique<Density>(placement, &densityGridCount_);

    printf("Initializing objective function system\n");
    printf("Placement region: (%.1f,%.1f) to (%.1f,%.1f)\n",
           placement.boundryLeft(),
           placement.boundryBottom(),
           placement.boundryRight(),
           placement.boundryTop());
}

double ObjectiveFunction::calculate() {
    const double wireValue = wirelength();
    const double densityValue = density();
    value_ = wireValue + lambda_ * densityValue;
    return value_;
}

void ObjectiveFunction::Backward() {
    wirelength_->Backward();
    density_->Backward();

#pragma omp parallel for schedule(static)
    for (size_t index = 0; index < grad_.size(); ++index) {
        grad_[index] = wirelength_->grad_[index] + lambda_ * density_->grad_[index];
    }
}

void ObjectiveFunction::initLambda() {
    wirelength_->Backward();
    density_->Backward();

    double wireMagnitude = 0.0;
    double densityMagnitude = 0.0;

#pragma omp parallel for reduction(+ : wireMagnitude, densityMagnitude) schedule(static)
    for (size_t index = 0; index < grad_.size(); ++index) {
        wireMagnitude += wirelength_->grad_[index].x * wirelength_->grad_[index].x +
                         wirelength_->grad_[index].y * wirelength_->grad_[index].y;
        densityMagnitude += density_->grad_[index].x * density_->grad_[index].x +
                            density_->grad_[index].y * density_->grad_[index].y;
    }

    if (densityMagnitude > kMinLambdaMagnitude) {
        lambda_ = std::sqrt(wireMagnitude) / std::sqrt(densityMagnitude);
    } else {
        lambda_ = 1.0;
    }
    printf("initial lambda %g\n", lambda_);
}

void ObjectiveFunction::updateLambda(double value) { lambda_ = value; }

void ObjectiveFunction::updateLambda() { lambda_ *= 2.0; }

void ObjectiveFunction::scaleLambda(double factor) { lambda_ *= factor; }

double ObjectiveFunction::wirelength() { return wirelength_->calculate(); }

double ObjectiveFunction::density() { return density_->calculate(); }

Placement& ObjectiveFunction::placement() { return *placement_; }
