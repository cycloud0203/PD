#include "floorplanner.h"
#include "sa_engine.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <climits>
#include <chrono>
#include <numeric>
#include <thread>

Floorplanner::Floorplanner(std::fstream& blockInput, std::fstream& netInput, double alpha)
    : _alpha(alpha), _outlineW(0), _outlineH(0), _numNets(0),
      _rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count())
{
    parseBlock(blockInput);
    parseNet(netInput);
    buildNetArrays();
    buildSharedData();
}

Floorplanner::~Floorplanner()
{
    for (auto b : _blocks) delete b;
    for (auto t : _terminals) delete t;
    for (auto n : _nets) delete n;
}

void Floorplanner::parseBlock(std::fstream& input)
{
    std::string token;
    int numBlocks = 0, numTerminals = 0;

    input >> token >> _outlineW >> _outlineH;
    input >> token >> numBlocks;
    input >> token >> numTerminals;

    for (int i = 0; i < numBlocks; ++i) {
        std::string name;
        int w, h;
        input >> name >> w >> h;
        _blocks.push_back(new Block(name, w, h));
        _nameToTerm[name] = _blocks.back();
    }

    for (int i = 0; i < numTerminals; ++i) {
        std::string name, keyword;
        int x, y;
        input >> name >> keyword >> x >> y;
        _terminals.push_back(new Terminal(name, x, y));
        _nameToTerm[name] = _terminals.back();
    }
}

void Floorplanner::parseNet(std::fstream& input)
{
    std::string token;
    int numNets = 0;
    input >> token >> numNets;

    for (int i = 0; i < numNets; ++i) {
        int degree;
        input >> token >> degree;
        Net* net = new Net();
        for (int j = 0; j < degree; ++j) {
            std::string name;
            input >> name;
            auto it = _nameToTerm.find(name);
            if (it != _nameToTerm.end())
                net->addTerm(it->second);
        }
        _nets.push_back(net);
    }
}

void Floorplanner::buildNetArrays()
{
    int nBlocks = (int)_blocks.size();

    std::map<std::string, int> nameToBlockId;
    for (int i = 0; i < nBlocks; ++i)
        nameToBlockId[_blocks[i]->getName()] = i;

    _numNets = (int)_nets.size();
    _netStart.resize(_numNets + 1);
    _allPins.clear();

    for (int i = 0; i < _numNets; ++i) {
        _netStart[i] = (int)_allPins.size();
        const auto& terms = _nets[i]->getTermList();
        for (const auto* term : terms) {
            NetPin pin;
            auto it = nameToBlockId.find(term->getName());
            if (it != nameToBlockId.end()) {
                pin.blockId = it->second;
                pin.fixedX2 = 0;
                pin.fixedY2 = 0;
            } else {
                pin.blockId = -1;
                pin.fixedX2 = term->getX1() + term->getX2();
                pin.fixedY2 = term->getY1() + term->getY2();
            }
            _allPins.push_back(pin);
        }
    }
    _netStart[_numNets] = (int)_allPins.size();

    _blockNets.resize(nBlocks);
    for (int i = 0; i < _numNets; ++i) {
        for (int j = _netStart[i]; j < _netStart[i + 1]; ++j) {
            if (_allPins[j].blockId >= 0)
                _blockNets[_allPins[j].blockId].push_back(i);
        }
    }
    for (auto& nets : _blockNets) {
        std::sort(nets.begin(), nets.end());
        nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
    }

    _blockAdj.resize(nBlocks);
    for (int i = 0; i < _numNets; ++i) {
        std::vector<int> blockIds;
        for (int j = _netStart[i]; j < _netStart[i + 1]; ++j) {
            if (_allPins[j].blockId >= 0)
                blockIds.push_back(_allPins[j].blockId);
        }
        for (size_t a = 0; a < blockIds.size(); ++a) {
            for (size_t b = a + 1; b < blockIds.size(); ++b) {
                _blockAdj[blockIds[a]].push_back(blockIds[b]);
                _blockAdj[blockIds[b]].push_back(blockIds[a]);
            }
        }
    }
    for (auto& adj : _blockAdj) {
        std::sort(adj.begin(), adj.end());
        std::vector<std::pair<int,int>> weighted;
        int idx = 0;
        while (idx < (int)adj.size()) {
            int j = idx;
            while (j < (int)adj.size() && adj[j] == adj[idx]) j++;
            weighted.push_back({-(j - idx), adj[idx]});
            idx = j;
        }
        std::sort(weighted.begin(), weighted.end());
        adj.clear();
        for (auto& [w, id] : weighted) adj.push_back(id);
    }
}

void Floorplanner::buildSharedData()
{
    int n = (int)_blocks.size();
    _sd.numBlocks = n;
    _sd.numNets   = _numNets;
    _sd.outlineW  = _outlineW;
    _sd.outlineH  = _outlineH;
    _sd.alpha     = _alpha;

    _sd.blockW.resize(n);
    _sd.blockH.resize(n);
    for (int i = 0; i < n; ++i) {
        _sd.blockW[i] = _blocks[i]->getWidth(false);
        _sd.blockH[i] = _blocks[i]->getHeight(false);
    }

    _sd.netStart  = std::move(_netStart);
    _sd.allPins   = std::move(_allPins);
    _sd.blockNets = std::move(_blockNets);
    _sd.blockAdj  = std::move(_blockAdj);

    _sd.normA = 1.0;
    _sd.normW = 1.0;
}

void Floorplanner::computeNormalization()
{
    int n = _sd.numBlocks;
    int m = std::max(5000, n * 100);
    double sumA = 0, sumW = 0;

    PerturbConfig cfg = PerturbConfig::forSize(n);
    BTree tempTree(n, _sd.blockW.data(), _sd.blockH.data(), _rng);
    WLCache wlc(n, _sd.numNets);
    tempTree.pack();

    for (int i = 0; i < m; ++i) {
        tempTree.perturb(cfg);
        tempTree.pack();
        sumA += (double)tempTree.getArea();
        sumW += wlc.compute(tempTree, _sd);
    }

    _sd.normA = sumA / m;
    _sd.normW = sumW / m;
    if (_sd.normA < 1.0) _sd.normA = 1.0;
    if (_sd.normW < 1.0) _sd.normW = 1.0;
}

double Floorplanner::computeWirelength()
{
    double wl = 0.0;
    for (const auto* net : _nets) wl += net->calcHPWL();
    return wl;
}

void Floorplanner::floorplan()
{
    auto startTime = std::chrono::steady_clock::now();
    int n = _sd.numBlocks;

    double timeLimitSec;
    if (n <= 12)      timeLimitSec = 8.0;
    else if (n <= 35) timeLimitSec = 100.0;
    else              timeLimitSec = 250.0;

    computeNormalization();

    unsigned hwThreads = std::thread::hardware_concurrency();
    if (hwThreads == 0) hwThreads = 4;
    int maxThreads;
    if (n <= 12)      maxThreads = 16;
    else if (n <= 35) maxThreads = 12;
    else              maxThreads = 8;
    int numThreads = std::min((int)hwThreads, maxThreads);

    auto deadline = startTime
                  + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                        std::chrono::duration<double>(timeLimitSec * 0.95));

    // Centralized multi-start sampling
    struct Candidate {
        double cost;
        std::vector<BTreeNode> nodes;
        int root;
    };
    int topK = std::min(numThreads, (n <= 12) ? 6 : (n <= 35) ? 4 : 3);
    std::vector<Candidate> topCandidates;

    {
        double sampleBudget = timeLimitSec * 0.03;
        BTree sampleTree(n, _sd.blockW.data(), _sd.blockH.data(), _rng);
        WLCache sampleWL(n, _sd.numNets);
        auto sampleEnd = std::chrono::steady_clock::now()
                       + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                             std::chrono::duration<double>(sampleBudget));
        int s = 0;
        while (std::chrono::steady_clock::now() < sampleEnd) {
            switch (s % 4) {
            case 0: sampleTree.init(); break;
            case 1: sampleTree.randomInit(); break;
            case 2: sampleTree.randomInsertionInit(); break;
            case 3: sampleTree.netOrderedInit(_sd.blockAdj); break;
            }
            sampleTree.pack();
            double area = (double)sampleTree.getArea();
            double wl = sampleWL.compute(sampleTree, _sd);
            double c = _sd.alpha * (area / _sd.normA)
                     + (1.0 - _sd.alpha) * (wl / _sd.normW);
            int exW = sampleTree.getWidth() - _sd.outlineW;
            int exH = sampleTree.getHeight() - _sd.outlineH;
            if (exW > 0 || exH > 0) {
                double rW = std::max(0.0, (double)exW / _sd.outlineW);
                double rH = std::max(0.0, (double)exH / _sd.outlineH);
                c *= (1.0 + 25.0 * (rW + rH));
            }

            bool insert = false;
            if ((int)topCandidates.size() < topK) {
                insert = true;
            } else if (c < topCandidates.back().cost) {
                insert = true;
            }
            if (insert) {
                topCandidates.push_back(
                    {c, sampleTree.getNodes(), sampleTree.getRoot()});
                std::sort(topCandidates.begin(), topCandidates.end(),
                    [](const Candidate& a, const Candidate& b) {
                        return a.cost < b.cost;
                    });
                if ((int)topCandidates.size() > topK)
                    topCandidates.resize(topK);
            }

            bool feas = (sampleTree.getWidth() <= _sd.outlineW
                      && sampleTree.getHeight() <= _sd.outlineH);
            double rc = _sd.alpha * (area / _sd.normA)
                      + (1.0 - _sd.alpha) * (wl / _sd.normW);
            _gb.tryUpdate(rc, feas,
                          sampleTree.getNodes(), sampleTree.getRoot(),
                          sampleTree.getWidth(), sampleTree.getHeight(),
                          sampleTree.blockXArr(), sampleTree.blockYArr(),
                          sampleTree.blockX2Arr(), sampleTree.blockY2Arr(), n);
            ++s;
        }
    }

    unsigned baseSeed = _rng();

    std::vector<std::thread> workers;
    for (int t = 0; t < numThreads; ++t) {
        if (t < (int)topCandidates.size()) {
            auto nodes = topCandidates[t].nodes;
            int root = topCandidates[t].root;
            workers.emplace_back([this, deadline, baseSeed, t,
                                  nodes = std::move(nodes), root]() {
                SAEngine engine(_sd, _gb, baseSeed + (unsigned)t, t);
                engine.runFromState(deadline, nodes, root);
            });
        } else {
            workers.emplace_back([this, deadline, baseSeed, t]() {
                SAEngine engine(_sd, _gb, baseSeed + (unsigned)t, t);
                engine.run(deadline);
            });
        }
    }

    for (auto& w : workers) w.join();

    // Steepest-descent local search polish (main thread)
    auto elapsed = [&]() -> double {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - startTime).count();
    };

    double lsTime = timeLimitSec * 0.97 - elapsed();
    if (lsTime > 0.02 && _gb.root >= 0) {
        BTree lsTree(n, _sd.blockW.data(), _sd.blockH.data(), _rng);
        WLCache lsWL(n, _sd.numNets);
        lsTree.loadState(_gb.nodes, _gb.root);
        lsTree.pack();

        auto lsCost = [&](BTree& tree) -> double {
            double area = (double)tree.getArea();
            double wl = lsWL.compute(tree, _sd);
            double c = _sd.alpha * (area / _sd.normA)
                     + (1.0 - _sd.alpha) * (wl / _sd.normW);
            int exW = tree.getWidth() - _sd.outlineW;
            int exH = tree.getHeight() - _sd.outlineH;
            if (exW > 0 || exH > 0) {
                double rW = std::max(0.0, (double)exW / _sd.outlineW);
                double rH = std::max(0.0, (double)exH / _sd.outlineH);
                c *= (1.0 + 25.0 * (rW + rH));
            }
            return c;
        };

        auto updateGB = [&](BTree& tree) {
            double area = (double)tree.getArea();
            double wl = lsWL.compute(tree, _sd);
            double rc = _sd.alpha * (area / _sd.normA)
                      + (1.0 - _sd.alpha) * (wl / _sd.normW);
            bool feas = (tree.getWidth() <= _sd.outlineW
                      && tree.getHeight() <= _sd.outlineH);
            _gb.tryUpdate(rc, feas,
                          tree.getNodes(), tree.getRoot(),
                          tree.getWidth(), tree.getHeight(),
                          tree.blockXArr(), tree.blockYArr(),
                          tree.blockX2Arr(), tree.blockY2Arr(), n);
        };

        auto findBestRotMask = [&](BTree& tree, double& bestCost) -> int {
            int bestMask = -1;
            for (int mask = 0; mask < (1 << n); ++mask) {
                for (int i = 0; i < n; ++i) {
                    int bid = tree.getNodes()[i].blockId;
                    tree.setRotation(i, ((mask >> bid) & 1) != 0);
                }
                tree.pack();
                if (tree.getWidth() > _sd.outlineW
                    || tree.getHeight() > _sd.outlineH) continue;
                double nc = lsCost(tree);
                if (nc < bestCost - 1e-12) {
                    bestCost = nc;
                    bestMask = mask;
                }
            }
            return bestMask;
        };

        auto applyRotMask = [&](BTree& tree, int mask) {
            for (int i = 0; i < n; ++i) {
                int bid = tree.getNodes()[i].blockId;
                tree.setRotation(i, ((mask >> bid) & 1) != 0);
            }
            tree.pack();
        };

        auto getRotMask = [&](const BTree& tree) -> int {
            int mask = 0;
            for (int i = 0; i < n; ++i) {
                if (tree.getNodes()[i].rotated)
                    mask |= (1 << tree.getNodes()[i].blockId);
            }
            return mask;
        };

        double curCost = lsCost(lsTree);

        if (n <= 12) {
            // Initial full rotation optimization
            {
                int origMask = getRotMask(lsTree);
                double bestCost = curCost;
                int bestMask = findBestRotMask(lsTree, bestCost);
                if (bestMask >= 0) {
                    applyRotMask(lsTree, bestMask);
                    curCost = bestCost;
                    updateGB(lsTree);
                } else {
                    applyRotMask(lsTree, origMask);
                }
            }

            // Joint swap + rotation enumeration
            bool improved = true;
            while (improved && elapsed() < timeLimitSec * 0.97) {
                improved = false;
                double bestSwapCost = curCost;
                int bestA = -1, bestB = -1, bestSwapMask = -1;
                int origMask = getRotMask(lsTree);

                for (int a = 0; a < n && elapsed() < timeLimitSec * 0.97; ++a) {
                    for (int b = a + 1; b < n; ++b) {
                        lsTree.swapBlockIds(a, b);
                        double trialCost = bestSwapCost;
                        int trialMask = findBestRotMask(lsTree, trialCost);
                        if (trialMask >= 0 && trialCost < bestSwapCost - 1e-12) {
                            bestSwapCost = trialCost;
                            bestA = a; bestB = b;
                            bestSwapMask = trialMask;
                        }
                        lsTree.swapBlockIds(a, b);
                    }
                }

                applyRotMask(lsTree, origMask);

                if (bestA >= 0) {
                    lsTree.swapBlockIds(bestA, bestB);
                    applyRotMask(lsTree, bestSwapMask);
                    curCost = bestSwapCost;
                    updateGB(lsTree);
                    improved = true;
                }
            }
        } else {
            bool improved = true;
            while (improved && elapsed() < timeLimitSec * 0.97) {
                improved = false;
                double bestDelta = 0;
                int bestA = -1, bestB = -1;
                bool bestIsRotation = false;

                for (int a = 0; a < n && elapsed() < timeLimitSec * 0.97; ++a) {
                    for (int b = a + 1; b < n; ++b) {
                        lsTree.swapBlockIds(a, b);
                        lsTree.pack();
                        double nc = lsCost(lsTree);
                        if (nc - curCost < bestDelta - 1e-12) {
                            bestDelta = nc - curCost;
                            bestA = a; bestB = b;
                            bestIsRotation = false;
                        }
                        lsTree.swapBlockIds(a, b);
                    }
                }

                for (int i = 0; i < n; ++i) {
                    lsTree.rotateNode(i);
                    lsTree.pack();
                    double nc = lsCost(lsTree);
                    if (nc - curCost < bestDelta - 1e-12) {
                        bestDelta = nc - curCost;
                        bestA = i;
                        bestIsRotation = true;
                    }
                    lsTree.rotateNode(i);
                }

                if (bestDelta < -1e-12) {
                    if (bestIsRotation)
                        lsTree.rotateNode(bestA);
                    else
                        lsTree.swapBlockIds(bestA, bestB);
                    lsTree.pack();
                    curCost = lsCost(lsTree);
                    updateGB(lsTree);
                    improved = true;
                }
            }
        }
    }
}

void Floorplanner::writeResult(std::fstream& output, double runtime)
{
    int n = (int)_blocks.size();

    if (_gb.bestX.empty() || _gb.root < 0) {
        output << 0 << std::endl;
        output << 0 << std::endl;
        output << 0 << std::endl;
        output << "0 0" << std::endl;
        output << std::fixed;
        output.precision(6);
        output << runtime << std::endl;
        for (int i = 0; i < n; ++i)
            output << _blocks[i]->getName() << " 0 0 0 0" << std::endl;
        return;
    }

    for (int i = 0; i < n; ++i)
        _blocks[i]->setPos(_gb.bestX[i], _gb.bestY[i],
                           _gb.bestX2[i], _gb.bestY2[i]);

    double totalWL = computeWirelength();
    long long area = (long long)_gb.width * _gb.height;
    double cost = _alpha * (double)area + (1.0 - _alpha) * totalWL;

    output << (long long)std::round(cost) << std::endl;
    output << (long long)std::round(totalWL) << std::endl;
    output << area << std::endl;
    output << _gb.width << " " << _gb.height << std::endl;
    output << std::fixed;
    output.precision(6);
    output << runtime << std::endl;

    for (int i = 0; i < n; ++i)
        output << _blocks[i]->getName() << " "
               << _gb.bestX[i] << " " << _gb.bestY[i] << " "
               << _gb.bestX2[i] << " " << _gb.bestY2[i] << std::endl;
}
