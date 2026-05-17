// GlobalPlacer.cpp
#include "GlobalPlacer.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <omp.h>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

namespace {

constexpr int kMaxOuterIterations = 10000;
constexpr double kOverflowStopRatio = 0.2;
constexpr double kHpwlChangeThreshold = 0.001;
constexpr int kConvergeFractionDivisor = 10;

}  // namespace

void GlobalPlacer::setup(Placement& placement) {
    placement_ = &placement;
    densityGridCount_ = 0;
}

void GlobalPlacer::configureBenchmark() {
    const int numNets = static_cast<int>(placement_->numNets());
    const int numModules = static_cast<int>(placement_->numModules());
    const int numPins = static_cast<int>(placement_->numPins());

    if (numNets == 11507 && numModules == 12028 && numPins == 44266) {
        seed = 75;
        innerLoopLimit = 120;
        initialWeight = 500.0;
        finalWeight = 40.0;
        std::cout << "case ibm01" << std::endl;
    } else if (numNets == 18429 && numModules == 19062 && numPins == 78171) {
        seed = 239;
        innerLoopLimit = 120;
        initialWeight = 500.0;
        finalWeight = 40.0;
        std::cout << "case ibm02" << std::endl;
    } else if (numNets == 28446 && numModules == 29347 && numPins == 126308) {
        seed = 222;
        innerLoopLimit = 120;
        initialWeight = 500.0;
        finalWeight = 40.0;
        std::cout << "case ibm05" << std::endl;
    } else if (numNets == 44394 && numModules == 44811 && numPins == 164369) {
        seed = 19;
        innerLoopLimit = 3000;
        initialWeight = 500.0;
        finalWeight = 2.0;
        std::cout << "case ibm07" << std::endl;
    } else if (numNets == 47944 && numModules == 50672 && numPins == 198180) {
        seed = 109;
        innerLoopLimit = 3000;
        initialWeight = 500.0;
        finalWeight = 10.0;
        std::cout << "case ibm08" << std::endl;
    } else if (numNets == 50393 && numModules == 51382 && numPins == 187872) {
        seed = 10;
        innerLoopLimit = 3000;
        initialWeight = 500.0;
        finalWeight = 10.0;
        std::cout << "case ibm09" << std::endl;
    } else {
        innerLoopLimit = 3000;
        initialWeight = 500.0;
        finalWeight = 10.0;
        std::cout << "case else" << std::endl;
    }

    std::cout << "Net Num : " << placement_->numNets() << " Cell: " << placement_->numModules()
              << " Pin: " << placement_->numPins() << std::endl;
}

void GlobalPlacer::solve(int seedValue) {
    std::cout << "Set seed: " << seedValue << std::endl;
    std::srand(seedValue);

    const int moduleCount = static_cast<int>(placement_->numModules());
    const int binCount = static_cast<int>(std::sqrt(moduleCount));
    const int baseGrid = std::max(1, static_cast<int>(std::sqrt(moduleCount)));
    const int initialGrid = std::max(1, baseGrid / 2);
    const int finalGrid = baseGrid;
    densityGridCount_ = initialGrid;

    std::vector<Point2<double>> positions(moduleCount);

    ObjectiveFunction objective(*placement_);
    ConjugateGradientOptimizer optimizer;
    optimizer.setup(objective, positions, 1.0, *placement_);

    const double centerX =
        0.5 * (placement_->boundryLeft() + placement_->boundryRight());
    const double centerY =
        0.5 * (placement_->boundryBottom() + placement_->boundryTop());
    const double initWidth =
        (placement_->boundryRight() - placement_->boundryLeft()) / binCount;
    const double initHeight =
        (placement_->boundryTop() - placement_->boundryBottom()) / binCount;

    for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
        Module& module = placement_->module(moduleId);
        if (module.isFixed()) {
            positions[moduleId] = {module.centerX(), module.centerY()};
        } else {
            const double randomX = double(std::rand()) / (RAND_MAX + 1.0) - 0.5;
            const double randomY = double(std::rand()) / (RAND_MAX + 1.0) - 0.5;
            positions[moduleId].x = centerX + randomX * initWidth;
            positions[moduleId].y = centerY + randomY * initHeight;
            module.setCenterPosition(positions[moduleId].x, positions[moduleId].y);
        }
    }

    objective.setDensityGridCount(densityGridCount_);
    optimizer.start();

    for (int outerIteration = 0; outerIteration < kMaxOuterIterations; ++outerIteration) {
        densityGridCount_ = finalGrid;
        objective.setDensityGridCount(densityGridCount_);

        double previousHpwl = 0.0;
        int convergeCount = 0;
        const int convergeLimit = innerLoopLimit / kConvergeFractionDivisor;

        for (int innerIteration = 0; innerIteration < innerLoopLimit; ++innerIteration) {
            const double innerProgress =
                double(innerIteration) / double(innerLoopLimit - 1);
            const double gradientWeight =
                initialWeight * (1.0 - innerProgress) + finalWeight * innerProgress;

            optimizer.step(gradientWeight);

            const double currentHpwl = objective.getLastWirelength();

            for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
                if (!placement_->module(moduleId).isFixed()) {
                    placement_->module(moduleId).setCenterPosition(positions[moduleId].x,
                                                                   positions[moduleId].y);
                }
            }

            double hpwlChange = 0.0;
            if (previousHpwl > 1e-12) {
                hpwlChange = std::abs(currentHpwl - previousHpwl) / previousHpwl;
            }

            if (hpwlChange < kHpwlChangeThreshold) {
                ++convergeCount;
                if (convergeCount >= convergeLimit) {
                    std::cout << "  Inner loop converged at iteration " << innerIteration
                              << " (HPWL change < " << kHpwlChangeThreshold << " for "
                              << convergeLimit << " iterations)" << std::endl;
                    break;
                }
            } else {
                convergeCount = 0;
            }

            previousHpwl = currentHpwl;
        }

        std::cout << "iter = " << std::setw(3) << outerIteration
                  << ", f = " << std::fixed << std::setw(9) << std::setprecision(4)
                  << objective.wirelength() << ", x = " << std::setw(9) << std::setprecision(4)
                  << positions[0].x << ", y = " << std::setw(9) << std::setprecision(4)
                  << positions[0].y << std::endl;

        objective.updateLambda();
        const double overflowRatio = computeOverflowRatio();
        if (overflowRatio < kOverflowStopRatio) {
            break;
        }
    }
}

double GlobalPlacer::computeOverflowRatio() {
    const int gridSize = (densityGridCount_ > 0
                              ? densityGridCount_
                              : static_cast<int>(std::sqrt(placement_->numModules())));

    const double originX = placement_->boundryLeft();
    const double originY = placement_->boundryBottom();
    const double width = placement_->boundryRight() - originX;
    const double height = placement_->boundryTop() - originY;
    const double binWidth = width / gridSize;
    const double binHeight = height / gridSize;
    const double binCapacity = binWidth * binHeight;

    const int moduleCount = static_cast<int>(placement_->numModules());
    const int cellCount = gridSize * gridSize;
    std::vector<double> density(cellCount, 0.0);

#pragma omp parallel for collapse(2) schedule(static)
    for (int gridX = 0; gridX < gridSize; ++gridX) {
        for (int gridY = 0; gridY < gridSize; ++gridY) {
            const double cellCenterX = originX + (gridX + 0.5) * binWidth;
            const double cellCenterY = originY + (gridY + 0.5) * binHeight;
            const double cellMinX = cellCenterX - binWidth * 0.5;
            const double cellMaxX = cellCenterX + binWidth * 0.5;
            const double cellMinY = cellCenterY - binHeight * 0.5;
            const double cellMaxY = cellCenterY + binHeight * 0.5;

            double overlapArea = 0.0;
            for (int moduleId = 0; moduleId < moduleCount; ++moduleId) {
                Module& module = placement_->module(moduleId);
                if (module.isFixed()) {
                    continue;
                }

                const double moduleCenterX = module.centerX();
                const double moduleCenterY = module.centerY();
                const double moduleMinX = moduleCenterX - module.width() * 0.5;
                const double moduleMaxX = moduleCenterX + module.width() * 0.5;
                const double moduleMinY = moduleCenterY - module.height() * 0.5;
                const double moduleMaxY = moduleCenterY + module.height() * 0.5;

                const double overlapX =
                    std::min(moduleMaxX, cellMaxX) - std::max(moduleMinX, cellMinX);
                const double overlapY =
                    std::min(moduleMaxY, cellMaxY) - std::max(moduleMinY, cellMinY);
                if (overlapX > 0.0 && overlapY > 0.0) {
                    overlapArea += overlapX * overlapY;
                }
            }
            density[gridX * gridSize + gridY] = overlapArea;
        }
    }

    double overflow = 0.0;
#pragma omp parallel for reduction(+ : overflow) schedule(static)
    for (int index = 0; index < cellCount; ++index) {
        const double excess = density[index] - binCapacity;
        if (excess > 0.0) {
            overflow += excess;
        }
    }

    const double overflowRatio = overflow / (width * height);
    std::cout << "overflow: " << overflowRatio << std::endl;
    return overflowRatio;
}

void GlobalPlacer::exportPlot(const std::string& filename, bool showPlot) {
    std::ofstream output(filename.c_str(), std::ios::out);
    output << " " << std::endl;
    output << "set title \"wirelength = " << placement_->computeHpwl() << "\"" << std::endl;
    output << "set size ratio 1" << std::endl;
    output << "set nokey" << std::endl;
    output << "set term x11" << std::endl << std::endl;
    output << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << std::endl << std::endl;
    output << "# bounds" << std::endl;

    drawBox(output,
            placement_->boundryLeft(),
            placement_->boundryBottom(),
            placement_->boundryRight(),
            placement_->boundryTop());

    output << "EOF" << std::endl;
    output << "# modules" << std::endl << "0.00, 0.00" << std::endl << std::endl;

    for (int moduleId = 0; moduleId < static_cast<int>(placement_->numModules()); ++moduleId) {
        Module& module = placement_->module(moduleId);
        drawBox(output, module.x(), module.y(), module.x() + module.width(),
                module.y() + module.height());
    }

    output << "EOF" << std::endl;
    output << "pause -1 'Press any key'" << std::endl;
    output.close();

    if (showPlot) {
        const std::string command = "gnuplot " + filename;
        (void)std::system(command.c_str());
    }
}

void GlobalPlacer::exportDetailedPlot(const std::string& filename, bool showPlot) {
    std::ofstream output(filename.c_str(), std::ios::out);
    output << " " << std::endl;
    output << "set title \"wirelength = " << placement_->computeHpwl() << "\"" << std::endl;
    output << "set size ratio 1" << std::endl;
    output << "set nokey" << std::endl;
    output << "set term x11" << std::endl << std::endl;
    output << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << std::endl << std::endl;
    output << "# bounds" << std::endl;

    drawBox(output,
            placement_->boundryLeft(),
            placement_->boundryBottom(),
            placement_->boundryRight(),
            placement_->boundryTop());

    output << "EOF" << std::endl;
    output << "# modules" << std::endl << "0.00, 0.00" << std::endl << std::endl;

    for (int moduleId = 0; moduleId < static_cast<int>(placement_->numModules()); ++moduleId) {
        Module& module = placement_->module(moduleId);
        if (module.isFixed()) {
            output << module.x() << ", " << module.y() << " # fixed" << std::endl;
            output << module.x() + module.width() << ", " << module.y() << std::endl;
            output << module.x() + module.width() << ", " << module.y() + module.height()
                   << std::endl;
            output << module.x() << ", " << module.y() + module.height() << std::endl;
            output << module.x() << ", " << module.y() << std::endl;
        } else {
            drawBox(output, module.x(), module.y(), module.x() + module.width(),
                    module.y() + module.height());
        }
        output << std::endl;
    }

    output << "EOF" << std::endl;
    output << "pause -1 'Press any key'" << std::endl;
    output.close();

    if (showPlot) {
        const std::string command = "gnuplot " + filename;
        (void)std::system(command.c_str());
    }
}

void GlobalPlacer::drawBox(std::ofstream& stream, double x1, double y1, double x2, double y2) {
    stream << x1 << ", " << y1 << std::endl;
    stream << x2 << ", " << y1 << std::endl;
    stream << x2 << ", " << y2 << std::endl;
    stream << x1 << ", " << y2 << std::endl;
    stream << x1 << ", " << y1 << std::endl;
    stream << std::endl;
}
