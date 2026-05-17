// Optimizer.h
#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <vector>

#include "ObjectiveFunction.h"
#include "Point.h"

class ConjugateGradientOptimizer {
public:
    void setup(ObjectiveFunction& objective,
               std::vector<Point2<double>>& variables,
               double initialStepSize,
               Placement& placement);

    void start();
    void step(double gradientWeight);

private:
    ObjectiveFunction* objective_;
    std::vector<Point2<double>>* variables_;
    Placement* placement_;

    std::vector<Point2<double>> prevGradient_;
    std::vector<Point2<double>> prevDirection_;
    int iterationCount_;
    double stepSize_;
};

#endif  // OPTIMIZER_H
