// ObjectiveFunction.h
#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <memory>
#include <vector>

#include "Placement.h"
#include "Point.h"

class BaseFunction {
public:
    std::vector<Point2<double>> grad_;
    double value_;

    explicit BaseFunction(size_t size) : value_(0.0) { grad_.resize(size); }
    virtual ~BaseFunction() = default;

    virtual double calculate() = 0;
    virtual void Backward() = 0;
};

class Wirelength : public BaseFunction {
public:
    explicit Wirelength(Placement& placement);

    double calculate() override;
    void Backward() override;

private:
    Placement* placement_;
    double gamma_;
};

class Density : public BaseFunction {
public:
    Density(Placement& placement, int* gridCountPtr);

    double calculate() override;
    void Backward() override;

private:
    Placement* placement_;
    int* gridCountPtr_;
};

class ObjectiveFunction : public BaseFunction {
public:
    explicit ObjectiveFunction(Placement& placement);
    ~ObjectiveFunction() override = default;

    double calculate() override;
    void Backward() override;

    double wirelength();
    double density();
    Placement& placement();

    double getLastWirelength() const { return wirelength_->value_; }

    void initLambda();
    void updateLambda(double value);
    void updateLambda();

    void setDensityGridCount(int count) { densityGridCount_ = count; }
    int densityGridCount() const { return densityGridCount_; }

private:
    std::unique_ptr<Wirelength> wirelength_;
    std::unique_ptr<Density> density_;
    Placement* placement_;
    int densityGridCount_;
    double lambda_;
};

#endif  // OBJECTIVEFUNCTION_H
