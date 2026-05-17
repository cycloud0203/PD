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
    // Per-benchmark density target factor for the GP penalty function.
    // Should sit at or below the legalizer's target density to leave a small
    // margin. Default 0.80 is safe for the 0.85 utilization cases.
    double densityTargetFactor = 0.80;

    // Coarse-grid overflow tracking. When > 0, computeCoarseOverflowRatio()
    // measures overflow at this coarser bin scale (sized to ~10 row heights
    // per side, matching check_density_target.pl). The solve() outer loop
    // uses it as a secondary early-exit criterion, which steers the optimizer
    // toward placements that legalize with less HPWL inflation.
    int coarseDensityGridCount = 0;

    void setup(Placement& placement);
    void configureBenchmark();
    void solve(int seedValue);
    double computeOverflowRatio();
    double computeCoarseOverflowRatio();

    void exportPlot(const std::string& filename, bool showPlot = false);
    void exportDetailedPlot(const std::string& filename, bool showPlot = false);

private:
    Placement* placement_;
    int densityGridCount_;

    void drawBox(std::ofstream& stream, double x1, double y1, double x2, double y2);
};

#endif  // GLOBALPLACER_H
