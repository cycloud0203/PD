#include "sa_engine.h"
#include "global_best.h"
#include <cmath>
#include <algorithm>

SAEngine::SAEngine(const SharedData& sd, GlobalBest& gb, unsigned seed, int initType)
    : _sd(sd), _gb(gb), _rng(seed),
      _tree(sd.numBlocks, sd.blockW.data(), sd.blockH.data(), _rng),
      _wlCache(sd.numBlocks, sd.numNets),
      _cfg(PerturbConfig::forSize(sd.numBlocks)),
      _initType(initType),
      _localBestCost(1e18), _localFeasible(false), _localBestRoot(-1)
{
    _tree.setBlockAdj(&sd.blockAdj);
}

double SAEngine::timeRemaining() const
{
    return std::chrono::duration<double>(
        _deadline - std::chrono::steady_clock::now()).count();
}

void SAEngine::run(std::chrono::steady_clock::time_point deadline)
{
    _deadline = deadline;
    if (timeRemaining() < 0.1) return;

    switch (_initType % 4) {
    case 0: _tree.init(); break;
    case 1: _tree.randomInit(); break;
    case 2: _tree.randomInsertionInit(); break;
    case 3: _tree.netOrderedInit(_sd.blockAdj); break;
    }
    _tree.pack();
    _wlCache.invalidate();
    double dA, dW;
    evalCostFull(dA, dW);
    _wlCache.accept();
    tryUpdateLocal(dA, dW);
    tryUpdateGlobal(dA, dW);

    runMultipleSA();
}

void SAEngine::runFromState(std::chrono::steady_clock::time_point deadline,
                            const std::vector<BTreeNode>& startNodes, int startRoot)
{
    _deadline = deadline;
    if (timeRemaining() < 0.1) return;

    _tree.loadState(startNodes, startRoot);
    _tree.pack();
    _wlCache.invalidate();
    double dA, dW;
    evalCostFull(dA, dW);
    _wlCache.accept();
    tryUpdateLocal(dA, dW);
    tryUpdateGlobal(dA, dW);

    runMultipleSA();
}

void SAEngine::runMultipleSA()
{
    int n = _sd.numBlocks;
    double totalBudget = timeRemaining();

    double wlRefBudget = 0;
    if (n <= 12) wlRefBudget = totalBudget * 0.25;

    double saBudget = totalBudget - wlRefBudget;

    double stage1Fraction = _localFeasible ? 0.0 : 0.40;
    double stage1Budget = saBudget * stage1Fraction;
    double stage2Budget = saBudget - stage1Budget;

    double perRunS1;
    if (n <= 12)      perRunS1 = stage1Budget * 0.40;
    else if (n <= 35) perRunS1 = stage1Budget * 0.30;
    else              perRunS1 = stage1Budget * 0.35;

    // Stage 1: Packing-Centric
    if (stage1Budget > 0.05 && !_localFeasible) {
        auto s1Deadline = std::chrono::steady_clock::now()
                        + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                              std::chrono::duration<double>(stage1Budget));

        int runIdx = 0;
        while (std::chrono::steady_clock::now() < s1Deadline) {
            double remaining = std::chrono::duration<double>(
                s1Deadline - std::chrono::steady_clock::now()).count();
            if (remaining < 0.05) break;

            double budget = std::min(perRunS1, remaining * 0.98);
            if (budget < 0.03) break;

            if (runIdx > 0) {
                bool fromBest = (runIdx % 3 == 2) && (_localBestRoot >= 0);
                if (fromBest) {
                    _tree.loadState(_localBestNodes, _localBestRoot);
                } else {
                    int initTypes = (runIdx + _initType) % 3;
                    if (initTypes == 0) _tree.randomInit();
                    else if (initTypes == 1) _tree.randomInsertionInit();
                    else _tree.netOrderedInit(_sd.blockAdj);
                }
                _tree.pack();
                _wlCache.invalidate();
                double dA, dW;
                evalCostFull(dA, dW);
                _wlCache.accept();
                tryUpdateLocal(dA, dW);
                tryUpdateGlobal(dA, dW);
            }

            bool feasible = runStageOne(budget);
            ++runIdx;
            if (feasible || _localFeasible) break;
        }
    }

    // Stage 2: Wirelength-Centric
    double perRunS2;
    if (n <= 12)      perRunS2 = stage2Budget * 0.28;
    else if (n <= 35) perRunS2 = stage2Budget * 0.20;
    else              perRunS2 = stage2Budget * 0.25;

    {
        double s2Remaining = std::min(stage2Budget,
            std::chrono::duration<double>(
                _deadline - std::chrono::steady_clock::now()).count() - wlRefBudget);
        if (s2Remaining < 0.05) s2Remaining = 0;

        auto s2Deadline = std::chrono::steady_clock::now()
                        + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                              std::chrono::duration<double>(s2Remaining));

        int runIdx = 0;
        while (std::chrono::steady_clock::now() < s2Deadline) {
            double remaining = std::chrono::duration<double>(
                s2Deadline - std::chrono::steady_clock::now()).count();
            if (remaining < 0.1) break;

            double budget = std::min(perRunS2, remaining * 0.98);
            if (budget < 0.05) break;

            if (runIdx > 0) {
                bool fromBest = (n <= 12) ? (runIdx % 2 == 1)
                                          : (runIdx % 3 == 2);
                if (fromBest && _localBestRoot >= 0) {
                    _tree.loadState(_localBestNodes, _localBestRoot);
                } else {
                    int initTypes = (runIdx + _initType) % 3;
                    if (initTypes == 0) _tree.randomInit();
                    else if (initTypes == 1) _tree.randomInsertionInit();
                    else _tree.netOrderedInit(_sd.blockAdj);
                }
                _tree.pack();
                _wlCache.invalidate();
                double dA, dW;
                evalCostFull(dA, dW);
                _wlCache.accept();
                tryUpdateLocal(dA, dW);
                tryUpdateGlobal(dA, dW);
            }

            runStageTwo(budget, _localFeasible);
            ++runIdx;
        }
    }

    if (wlRefBudget > 0.1 && _localBestRoot >= 0 && _localFeasible) {
        _tree.loadState(_localBestNodes, _localBestRoot);
        _tree.pack();
        _wlCache.invalidate();
        runWLRefinement(std::min(wlRefBudget, timeRemaining() * 0.95));
    }
}

double SAEngine::evalCostPacking(double& outArea)
{
    outArea = (double)_tree.getArea();
    double c = outArea / _sd.normA;

    int exW = _tree.getWidth() - _sd.outlineW;
    int exH = _tree.getHeight() - _sd.outlineH;
    if (exW > 0 || exH > 0) {
        double rW = std::max(0.0, (double)exW / _sd.outlineW);
        double rH = std::max(0.0, (double)exH / _sd.outlineH);
        c *= (1.0 + 50.0 * (rW + rH));
    }

    double targetAR = (double)_sd.outlineW / _sd.outlineH;
    double curAR = (_tree.getHeight() > 0)
                 ? (double)_tree.getWidth() / _tree.getHeight() : 1e9;
    double arDiff = std::abs(curAR - targetAR) / targetAR;
    c *= (1.0 + 2.0 * arDiff);

    return c;
}

double SAEngine::evalCostFull(double& outArea, double& outWL)
{
    outArea = (double)_tree.getArea();
    outWL = _wlCache.compute(_tree, _sd);
    double c = _sd.alpha * (outArea / _sd.normA)
             + (1.0 - _sd.alpha) * (outWL / _sd.normW);

    int exW = _tree.getWidth() - _sd.outlineW;
    int exH = _tree.getHeight() - _sd.outlineH;
    if (exW > 0 || exH > 0) {
        double rW = std::max(0.0, (double)exW / _sd.outlineW);
        double rH = std::max(0.0, (double)exH / _sd.outlineH);
        c *= (1.0 + 25.0 * (rW + rH));
    }

    double targetAR = (double)_sd.outlineW / _sd.outlineH;
    double curAR = (_tree.getHeight() > 0)
                 ? (double)_tree.getWidth() / _tree.getHeight() : 1e9;
    double arDiff = std::abs(curAR - targetAR) / targetAR;
    c *= (1.0 + 2.0 * arDiff);

    return c;
}

double SAEngine::computeInitTemp(bool packingOnly, int warmupSteps)
{
    _tree.pack();
    _wlCache.invalidate();
    double prevCost;
    if (packingOnly) {
        double dA;
        prevCost = evalCostPacking(dA);
    } else {
        double dA, dW;
        prevCost = evalCostFull(dA, dW);
    }

    double sumAbsDelta = 0;
    for (int i = 0; i < warmupSteps; ++i) {
        _tree.perturb(_cfg);
        _tree.pack();
        double newCost;
        if (packingOnly) {
            double dA;
            newCost = evalCostPacking(dA);
        } else {
            double dA, dW;
            newCost = evalCostFull(dA, dW);
        }
        sumAbsDelta += std::abs(newCost - prevCost);
        _tree.undoPerturb();
    }

    _wlCache.invalidate();

    double avgDelta = sumAbsDelta / warmupSteps;
    if (avgDelta < 1e-10) avgDelta = 1e-10;
    return -avgDelta / std::log(0.95);
}

void SAEngine::tryUpdateLocal(double area, double wl)
{
    double realCost = _sd.alpha * (area / _sd.normA)
                    + (1.0 - _sd.alpha) * (wl / _sd.normW);
    bool feas = (_tree.getWidth() <= _sd.outlineW
              && _tree.getHeight() <= _sd.outlineH);

    bool doUpdate = false;
    if (_localBestRoot == -1)                        doUpdate = true;
    else if (feas && !_localFeasible)                doUpdate = true;
    else if (feas == _localFeasible && realCost < _localBestCost) doUpdate = true;

    if (doUpdate) {
        _localBestCost = realCost;
        _localFeasible = feas;
        _localBestNodes = _tree.getNodes();
        _localBestRoot = _tree.getRoot();
    }
}

void SAEngine::tryUpdateGlobal(double area, double wl)
{
    double realCost = _sd.alpha * (area / _sd.normA)
                    + (1.0 - _sd.alpha) * (wl / _sd.normW);
    bool feas = (_tree.getWidth() <= _sd.outlineW
              && _tree.getHeight() <= _sd.outlineH);

    _gb.tryUpdate(realCost, feas,
                  _tree.getNodes(), _tree.getRoot(),
                  _tree.getWidth(), _tree.getHeight(),
                  _tree.blockXArr(), _tree.blockYArr(),
                  _tree.blockX2Arr(), _tree.blockY2Arr(),
                  _sd.numBlocks);
}

bool SAEngine::runStageOne(double budget)
{
    int n = _sd.numBlocks;

    double T = computeInitTemp(true, 300);
    _tree.pack();
    _wlCache.invalidate();

    double dA;
    double curCost = evalCostPacking(dA);

    int movesPerTemp;
    if (n <= 12) movesPerTemp = n * n * 8;
    else         movesPerTemp = std::max(n * n * 3, 80);

    double frozenTemp = T * 1e-7;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    auto saEnd = std::chrono::steady_clock::now()
               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(budget * 0.95));
    if (saEnd > _deadline) saEnd = _deadline;

    while (T > frozenTemp) {
        if (std::chrono::steady_clock::now() >= saEnd) break;

        int accepted = 0;
        for (int i = 0; i < movesPerTemp; ++i) {
            _tree.perturb(_cfg);
            _tree.pack();
            double newA;
            double newCost = evalCostPacking(newA);
            double delta = newCost - curCost;

            if (delta <= 0 || unif(_rng) < std::exp(-delta / T)) {
                curCost = newCost;
                ++accepted;

                bool feas = (_tree.getWidth() <= _sd.outlineW
                          && _tree.getHeight() <= _sd.outlineH);
                if (feas) {
                    _wlCache.invalidate();
                    double wl = _wlCache.compute(_tree, _sd);
                    _wlCache.accept();
                    tryUpdateLocal(newA, wl);
                    tryUpdateGlobal(newA, wl);
                    return true;
                }
            } else {
                _tree.undoPerturb();
            }
        }

        double acceptRate = (double)accepted / movesPerTemp;
        if (acceptRate > 0.8)       T *= 0.8;
        else if (acceptRate >= 0.15) T *= 0.995;
        else                         T *= 0.9;
    }

    return false;
}

void SAEngine::runStageTwo(double budget, bool hardConstraint)
{
    int n = _sd.numBlocks;

    double T = computeInitTemp(false, 500);
    _tree.pack();
    _wlCache.invalidate();

    double dA, dW;
    double curCost = evalCostFull(dA, dW);
    _wlCache.accept();

    int movesPerTemp;
    if (n <= 12) movesPerTemp = n * n * 8;
    else         movesPerTemp = std::max(n * n * 3, 80);

    double frozenTemp = T * 1e-7;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    double saTimeBudget = budget * 0.90;
    auto saEnd = std::chrono::steady_clock::now()
               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(saTimeBudget));
    if (saEnd > _deadline) saEnd = _deadline;

    while (T > frozenTemp) {
        if (std::chrono::steady_clock::now() >= saEnd) break;

        int accepted = 0;
        for (int i = 0; i < movesPerTemp; ++i) {
            _tree.perturb(_cfg);
            _tree.pack();

            if (hardConstraint &&
                (_tree.getWidth() > _sd.outlineW || _tree.getHeight() > _sd.outlineH)) {
                _tree.undoPerturb();
                continue;
            }

            double newA, newW;
            double newCost = evalCostFull(newA, newW);
            double delta = newCost - curCost;

            if (delta <= 0 || unif(_rng) < std::exp(-delta / T)) {
                curCost = newCost;
                _wlCache.accept();
                ++accepted;
                tryUpdateLocal(newA, newW);
                tryUpdateGlobal(newA, newW);
            } else {
                _tree.undoPerturb();
            }
        }

        double acceptRate = (double)accepted / movesPerTemp;
        if (acceptRate > 0.8)       T *= 0.8;
        else if (acceptRate >= 0.15) T *= 0.995;
        else                         T *= 0.9;
    }

    auto greedyEnd = std::chrono::steady_clock::now()
                   + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                         std::chrono::duration<double>(budget * 0.99 - saTimeBudget));
    if (greedyEnd > _deadline) greedyEnd = _deadline;

    while (std::chrono::steady_clock::now() < greedyEnd) {
        _tree.perturb(_cfg);
        _tree.pack();

        if (hardConstraint &&
            (_tree.getWidth() > _sd.outlineW || _tree.getHeight() > _sd.outlineH)) {
            _tree.undoPerturb();
            continue;
        }

        double gA, gW;
        double nc = evalCostFull(gA, gW);
        if (nc <= curCost) {
            curCost = nc;
            _wlCache.accept();
            tryUpdateLocal(gA, gW);
            tryUpdateGlobal(gA, gW);
        } else {
            _tree.undoPerturb();
        }
    }
}

void SAEngine::runSingleSA(double budget)
{
    int n = _sd.numBlocks;

    double T = computeInitTemp(false, 500);
    _tree.pack();
    _wlCache.invalidate();

    double dA, dW;
    double curCost = evalCostFull(dA, dW);
    _wlCache.accept();

    int movesPerTemp;
    if (n <= 12) movesPerTemp = n * n * 8;
    else         movesPerTemp = std::max(n * n * 3, 80);

    double frozenTemp = T * 1e-7;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    double saTimeBudget = budget * 0.90;
    auto saEnd = std::chrono::steady_clock::now()
               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(saTimeBudget));

    while (T > frozenTemp) {
        if (std::chrono::steady_clock::now() >= saEnd) break;

        int accepted = 0;
        for (int i = 0; i < movesPerTemp; ++i) {
            _tree.perturb(_cfg);
            _tree.pack();
            double newA, newW;
            double newCost = evalCostFull(newA, newW);
            double delta = newCost - curCost;

            if (delta <= 0 || unif(_rng) < std::exp(-delta / T)) {
                curCost = newCost;
                _wlCache.accept();
                ++accepted;
                tryUpdateLocal(newA, newW);
                tryUpdateGlobal(newA, newW);
            } else {
                _tree.undoPerturb();
            }
        }

        double acceptRate = (double)accepted / movesPerTemp;
        if (acceptRate > 0.8)       T *= 0.8;
        else if (acceptRate >= 0.15) T *= 0.995;
        else                         T *= 0.9;
    }

    auto greedyEnd = std::chrono::steady_clock::now()
                   + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                         std::chrono::duration<double>(budget * 0.99 - saTimeBudget));
    if (greedyEnd > _deadline) greedyEnd = _deadline;

    while (std::chrono::steady_clock::now() < greedyEnd) {
        _tree.perturb(_cfg);
        _tree.pack();
        double gA, gW;
        double nc = evalCostFull(gA, gW);
        if (nc <= curCost) {
            curCost = nc;
            _wlCache.accept();
            tryUpdateLocal(gA, gW);
            tryUpdateGlobal(gA, gW);
        } else {
            _tree.undoPerturb();
        }
    }
}

void SAEngine::runWLRefinement(double budget)
{
    int n = _sd.numBlocks;
    double refAlpha = _sd.alpha * 0.4;

    double T = computeInitTemp(false, 200) * 0.05;
    _tree.pack();
    _wlCache.invalidate();

    double wl = _wlCache.compute(_tree, _sd);
    _wlCache.accept();
    double area = (double)_tree.getArea();
    double curCost = refAlpha * (area / _sd.normA) + (1.0 - refAlpha) * (wl / _sd.normW);

    int movesPerTemp = std::max(n * n * 5, 100);

    auto saEnd = std::chrono::steady_clock::now()
               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(budget * 0.95));
    if (saEnd > _deadline) saEnd = _deadline;

    double frozenTemp = T * 1e-5;

    std::uniform_real_distribution<double> unif(0.0, 1.0);

    while (T > frozenTemp && std::chrono::steady_clock::now() < saEnd) {
        int accepted = 0;
        for (int i = 0; i < movesPerTemp; ++i) {
            _tree.perturb(_cfg);
            _tree.pack();

            if (_tree.getWidth() > _sd.outlineW || _tree.getHeight() > _sd.outlineH) {
                _tree.undoPerturb();
                continue;
            }

            double newWL = _wlCache.compute(_tree, _sd);
            double newArea = (double)_tree.getArea();
            double newCost = refAlpha * (newArea / _sd.normA)
                           + (1.0 - refAlpha) * (newWL / _sd.normW);
            double delta = newCost - curCost;

            if (delta <= 0 || unif(_rng) < std::exp(-delta / T)) {
                curCost = newCost;
                wl = newWL;
                area = newArea;
                _wlCache.accept();
                ++accepted;
                tryUpdateLocal(newArea, newWL);
                tryUpdateGlobal(newArea, newWL);
            } else {
                _tree.undoPerturb();
            }
        }

        double acceptRate = (double)accepted / movesPerTemp;
        if (acceptRate > 0.8)       T *= 0.8;
        else if (acceptRate >= 0.15) T *= 0.995;
        else                         T *= 0.9;
    }

    while (std::chrono::steady_clock::now() < saEnd) {
        _tree.perturb(_cfg);
        _tree.pack();
        if (_tree.getWidth() > _sd.outlineW || _tree.getHeight() > _sd.outlineH) {
            _tree.undoPerturb();
            continue;
        }
        double newWL = _wlCache.compute(_tree, _sd);
        double newArea = (double)_tree.getArea();
        double newCost = refAlpha * (newArea / _sd.normA)
                       + (1.0 - refAlpha) * (newWL / _sd.normW);
        if (newCost <= curCost) {
            curCost = newCost;
            wl = newWL;
            area = newArea;
            _wlCache.accept();
            tryUpdateLocal(newArea, newWL);
            tryUpdateGlobal(newArea, newWL);
        } else {
            _tree.undoPerturb();
        }
    }
}
