// Optimizer.cpp
#include "Optimizer.h"

#include <cmath>

void ConjugateGradientOptimizer::setup(ObjectiveFunction& objective,
                                       std::vector<Point2<double>>& variables,
                                       double initialStepSize,
                                       Placement& placement) {
    objective_ = &objective;
    variables_ = &variables;
    stepSize_ = initialStepSize;
    placement_ = &placement;
    iterationCount_ = 0;

    const size_t variableCount = variables.size();
    prevGradient_.resize(variableCount);
    prevDirection_.resize(variableCount);

    printf("Optimizer ready\n");
    printf("Area: (%.1f,%.1f) to (%.1f,%.1f)\n",
           placement.boundryLeft(),
           placement.boundryBottom(),
           placement.boundryRight(),
           placement.boundryTop());
}

void ConjugateGradientOptimizer::start() {
    iterationCount_ = 0;
    objective_->initLambda();
}

void ConjugateGradientOptimizer::step(double gradientWeight) {
    const int variableCount = static_cast<int>(variables_->size());
    objective_->Backward();

    double beta = 0.0;
    std::vector<Point2<double>> direction(variableCount);

    if (iterationCount_ == 0) {
        for (int index = 0; index < variableCount; ++index) {
            direction[index] = Point2<double>(
                -gradientWeight * objective_->grad_[index].x,
                -gradientWeight * objective_->grad_[index].y);
        }
    } else {
        double numerator = 0.0;
        double denominator = 0.0;
        for (int index = 0; index < variableCount; ++index) {
            const Point2<double>& gradient = objective_->grad_[index];
            const Point2<double>& previousGradient = prevGradient_[index];
            const Point2<double> gradientDelta(
                gradient.x - previousGradient.x, gradient.y - previousGradient.y);

            numerator += gradient.x * gradientDelta.x + gradient.y * gradientDelta.y;
            denominator += std::abs(gradient.x) + std::abs(gradient.y);
        }

        beta = numerator / (denominator * denominator);

        for (int index = 0; index < variableCount; ++index) {
            direction[index].x =
                -objective_->grad_[index].x * gradientWeight + beta * prevDirection_[index].x;
            direction[index].y =
                -objective_->grad_[index].y * gradientWeight + beta * prevDirection_[index].y;
        }
    }

    double directionNormSquared = 0.0;
    const int binCount = static_cast<int>(std::sqrt(placement_->numModules()));
    const double binWidth =
        (placement_->boundryRight() - placement_->boundryLeft()) / binCount;

    for (int index = 0; index < variableCount; ++index) {
        if (placement_->module(index).isFixed()) {
            continue;
        }
        directionNormSquared += direction[index].x * direction[index].x +
                                direction[index].y * direction[index].y;
    }

    const double directionNorm = std::sqrt(directionNormSquared);
    constexpr double kMinDirectionNorm = 1e-12;
    if (directionNorm > kMinDirectionNorm) {
        stepSize_ = 100.0 * binWidth / directionNorm;
    } else {
        stepSize_ = 0.0;
    }

    for (int index = 0; index < variableCount; ++index) {
        if (placement_->module(index).isFixed()) {
            continue;
        }

        (*variables_)[index].x += stepSize_ * direction[index].x;
        (*variables_)[index].y += stepSize_ * direction[index].y;

        const double halfWidth = placement_->module(index).width() * 0.5;
        const double halfHeight = placement_->module(index).height() * 0.5;
        const double maxX = objective_->placement().boundryRight() - halfWidth;
        const double maxY = objective_->placement().boundryTop() - halfHeight;
        const double minX = objective_->placement().boundryLeft() + halfWidth;
        const double minY = objective_->placement().boundryBottom() + halfHeight;

        if ((*variables_)[index].x > maxX) {
            (*variables_)[index].x = maxX;
        }
        if ((*variables_)[index].x < minX) {
            (*variables_)[index].x = minX;
        }
        if ((*variables_)[index].y > maxY) {
            (*variables_)[index].y = maxY;
        }
        if ((*variables_)[index].y < minY) {
            (*variables_)[index].y = minY;
        }
    }

    prevGradient_ = objective_->grad_;
    prevDirection_ = direction;
    ++iterationCount_;
}
