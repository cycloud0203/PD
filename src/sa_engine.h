#ifndef SA_ENGINE_H
#define SA_ENGINE_H

#include <chrono>
#include <random>
#include <vector>
#include "shared_data.h"
#include "btree.h"
#include "wl_cache.h"

struct GlobalBest;

class SAEngine {
public:
    SAEngine(const SharedData& sd, GlobalBest& gb, unsigned seed, int initType);
    void run(std::chrono::steady_clock::time_point deadline);
    void runFromState(std::chrono::steady_clock::time_point deadline,
                      const std::vector<BTreeNode>& startNodes, int startRoot);

private:
    const SharedData& _sd;
    GlobalBest& _gb;
    std::mt19937 _rng;
    BTree _tree;
    WLCache _wlCache;
    PerturbConfig _cfg;
    int _initType;
    std::chrono::steady_clock::time_point _deadline;

    double _localBestCost;
    bool _localFeasible;
    std::vector<BTreeNode> _localBestNodes;
    int _localBestRoot;

    double timeRemaining() const;

    double evalCostPacking(double& outArea);
    double evalCostFull(double& outArea, double& outWL);

    double computeInitTemp(bool packingOnly, int warmupSteps = 1000);

    void runMultipleSA();
    void runSingleSA(double budget);
    void runWLRefinement(double budget);
    bool runStageOne(double budget);
    void runStageTwo(double budget, bool hardConstraint);

    void tryUpdateLocal(double area, double wl);
    void tryUpdateGlobal(double area, double wl);
};

#endif
