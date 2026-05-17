#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include <fstream>
#include <string>

#include "Placement.h"

class GlobalPlacer {
public:
    int seed = 0;
    int innerLoopLimit = 3000;
    double initialWeight = 500.0;
    double finalWeight = 10.0;

    void setup(Placement& placement);
    void configureBenchmark();
    void solve(int seedValue);
    double computeOverflowRatio();

    void exportPlot(const std::string& filename, bool showPlot = false);
    void exportDetailedPlot(const std::string& filename, bool showPlot = false);

private:
    Placement* placement_;
    int densityGridCount_;

    void drawBox(std::ofstream& stream, double x1, double y1, double x2, double y2);
};

#endif  // GLOBALPLACER_H
