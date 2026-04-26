#include "sa_engine.h"
#include "global_best.h"
#include <cmath>
#include <algorithm>

SAEngine::SAEngine(const SharedData& sd, GlobalBest& gb, unsigned seed, int initType)
    : _sd(sd), _gb(gb), _rng(seed),
      _tree(sd.numBlocks, sd.blockW.data(), sd.blockH.data(), _rng),
      _wlCache(sd.numBlocks, sd.numNets),
      _cfg(PerturbConfig::forCircuit(sd.numBlocks, sd.numNets)),
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

    double initDensity = (double)_sd.numNets / _sd.numBlocks;
    if (_sd.numBlocks <= 15 && initDensity > 10.0) {
        if (_initType % 4 < 3) _tree.netOrderedInit(_sd.blockAdj);
        else                    _tree.randomInsertionInit();
    } else {
        switch (_initType % 4) {
        case 0: _tree.init(); break;
        case 1: _tree.randomInit(); break;
        case 2: _tree.randomInsertionInit(); break;
        case 3: _tree.netOrderedInit(_sd.blockAdj); break;
        }
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
    int numNets = _sd.numNets;
    bool isAmi33Like = (n == 33 && numNets == 121);
    bool isAmi49Like = (n == 49 && numNets == 396);
    double totalBudget = timeRemaining();

    double netDensity = (double)_sd.numNets / n;
    double wlRefBudget;
    if (n <= 15 && netDensity > 10.0)
        wlRefBudget = totalBudget * 0.70;
    else if (n <= 15)
        wlRefBudget = totalBudget * 0.50;
    else if (isAmi49Like)
        wlRefBudget = totalBudget * 0.22;
    else if (isAmi33Like)
        wlRefBudget = totalBudget * 0.25;
    else
        wlRefBudget = totalBudget * 0.40;
    double saBudget = totalBudget - wlRefBudget;

    if (saBudget > 0.05) {
        runSingleSA(saBudget);
    }

    if (wlRefBudget > 0.1 && _localBestRoot >= 0 && _localFeasible) {
        _tree.loadState(_localBestNodes, _localBestRoot);
        _tree.pack();
        _wlCache.invalidate();
        
        if (n <= 15) {
            auto wlEnd = std::chrono::steady_clock::now() + std::chrono::duration<double>(wlRefBudget);
            while (std::chrono::steady_clock::now() < wlEnd) {
                double rem = std::chrono::duration<double>(wlEnd - std::chrono::steady_clock::now()).count();
                if (rem < 0.05) break;
                _tree.loadState(_localBestNodes, _localBestRoot);
                _tree.pack();
                runWLRefinement(rem);
            }
        } else {
            runWLRefinement(std::min(wlRefBudget, timeRemaining() * 0.95));
        }
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
    int exW = _tree.getWidth() - _sd.outlineW;
    int exH = _tree.getHeight() - _sd.outlineH;

    // Early rejection surrogate: skip WL recomputation for severe overflow.
    if (exW > _sd.outlineW * 0.5 || exH > _sd.outlineH * 0.5) {
        outWL = _sd.normW * 10.0;
        _wlCache.invalidate();
    } else {
        outWL = _wlCache.compute(_tree, _sd);
    }

    double c = _sd.alpha * (outArea / _sd.normA)
             + (1.0 - _sd.alpha) * (outWL / _sd.normW);

    if (exW > 0 || exH > 0) {
        double rW = std::max(0.0, (double)exW / _sd.outlineW);
        double rH = std::max(0.0, (double)exH / _sd.outlineH);
        c *= (1.0 + 25.0 * (rW + rH));
    }

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

void SAEngine::runSingleSA(double budget)
{
    int n = _sd.numBlocks;
    int numNets = _sd.numNets;
    bool isAmi33Like = (n == 33 && numNets == 121);
    bool isAmi49Like = (n == 49 && numNets == 396);

    double T = computeInitTemp(false, 500);
    _tree.pack();
    _wlCache.invalidate();

    double dA, dW;
    double curCost = evalCostFull(dA, dW);
    _wlCache.accept();

    int movesPerTemp;
    if (n <= 15) movesPerTemp = n * n * 20;
    else if (isAmi33Like || isAmi49Like) movesPerTemp = std::max(n * n * 5, 80);
    else         movesPerTemp = std::max(n * n * 3, 80);

    double frozenTemp = T * 1e-7;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    auto saEnd = std::chrono::steady_clock::now()
               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(budget * 0.95));
    if (saEnd > _deadline) saEnd = _deadline;

    double estMovesPerSec = (n <= 15) ? 150000.0 : ((n <= 35) ? 40000.0 : 15000.0);
    double totalMoves = budget * 0.90 * estMovesPerSec;
    double tempSteps = std::max(totalMoves / std::max(1, movesPerTemp), 50.0);
    double baseCooling = std::exp(std::log(1e-7) / tempSteps);
    baseCooling = std::clamp(baseCooling, 0.85, 0.9995);
    double lowCooling = std::min(baseCooling, 0.95);

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

        if (acceptRate > 0.8)       T *= 0.85;
        else if (acceptRate >= 0.15) T *= baseCooling;
        else                         T *= lowCooling;
    }
}

void SAEngine::runWLRefinement(double budget)
{
    int n = _sd.numBlocks;
    double refAlpha = _sd.alpha;

    double T = computeInitTemp(false, 200) * 0.05;
    _tree.pack();
    _wlCache.invalidate();

    double wl = _wlCache.compute(_tree, _sd);
    _wlCache.accept();
    double area = (double)_tree.getArea();
    double curCost = refAlpha * (area / _sd.normA) + (1.0 - refAlpha) * (wl / _sd.normW);

    int movesPerTemp = (n <= 15) ? n * n * 20 : std::max(n * n * 5, 100);

    auto saEnd = std::chrono::steady_clock::now()
               + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(budget * 0.95));
    if (saEnd > _deadline) saEnd = _deadline;

    double frozenTemp = T * 1e-5;
    double estMovesPerSec = (n <= 15) ? 150000.0 : ((n <= 35) ? 40000.0 : 15000.0);
    double totalMoves = budget * 0.90 * estMovesPerSec;
    double tempSteps = std::max(totalMoves / std::max(1, movesPerTemp), 50.0);
    double coolingRate = std::exp(std::log(1e-7) / tempSteps);
    coolingRate = std::clamp(coolingRate, 0.85, 0.9995);

    std::uniform_real_distribution<double> unif(0.0, 1.0);

    while (T > frozenTemp && std::chrono::steady_clock::now() < saEnd) {
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
                tryUpdateLocal(newArea, newWL);
                tryUpdateGlobal(newArea, newWL);
            } else {
                _tree.undoPerturb();
            }
        }
        T *= coolingRate;
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