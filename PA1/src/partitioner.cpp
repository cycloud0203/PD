#include "partitioner.h"
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <chrono>
#include <random>
#include <queue>
#include <iostream>

using namespace std;

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

void Partitioner::parseInput(ifstream& inFile)
{
    string str;
    inFile >> str;
    _bFactor = stod(str);

    unordered_map<string, int> cellMap;
    cellMap.reserve(1 << 18);

    vector<int> netCellBuf;

    while (inFile >> str) {
        if (str != "NET") continue;

        string netName;
        inFile >> netName;
        int netId = _netNum++;
        _netCells.emplace_back();
        netCellBuf.clear();

        string tok;
        while (inFile >> tok) {
            if (tok == ";") break;

            auto it = cellMap.find(tok);
            int cellId;
            if (it == cellMap.end()) {
                cellId = _cellNum++;
                cellMap.emplace(tok, cellId);
                _cellName.push_back(tok);
                _cellNets.emplace_back();
            } else {
                cellId = it->second;
            }

            bool dup = false;
            for (int s : netCellBuf)
                if (s == cellId) { dup = true; break; }

            if (!dup) {
                netCellBuf.push_back(cellId);
                _cellNets[cellId].push_back(netId);
                _netCells[netId].push_back(cellId);
            }
        }
    }

    // [Opt 3] Remove single-cell nets: they can never be cut
    for (int n = 0; n < _netNum; ++n) {
        if ((int)_netCells[n].size() <= 1) {
            if (_netCells[n].size() == 1) {
                int c = _netCells[n][0];
                auto& cn = _cellNets[c];
                cn.erase(remove(cn.begin(), cn.end(), n), cn.end());
            }
            _netCells[n].clear();
        }
    }

    // [Opt 8] Sort adjacency lists by ID for better cache locality
    for (int n = 0; n < _netNum; ++n)
        sort(_netCells[n].begin(), _netCells[n].end());
    for (int c = 0; c < _cellNum; ++c)
        sort(_cellNets[c].begin(), _cellNets[c].end());

    _maxPinNum = 0;
    for (int i = 0; i < _cellNum; ++i)
        _maxPinNum = max(_maxPinNum, (int)_cellNets[i].size());

    _cells.resize(_cellNum);
    _nets.resize(_netNum);
}

// ---------------------------------------------------------------------------
// Initial partition helpers
// ---------------------------------------------------------------------------

void Partitioner::initPartition()
{
    int half = _cellNum / 2;
    _partSize[0] = half;
    _partSize[1] = _cellNum - half;
    for (int i = 0; i < _cellNum; ++i)
        _cells[i].part = (i < half) ? 0 : 1;
}

void Partitioner::randomPartition(mt19937& rng)
{
    vector<int> perm(_cellNum);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng);

    int half = _cellNum / 2;
    _partSize[0] = half;
    _partSize[1] = _cellNum - half;
    for (int i = 0; i < _cellNum; ++i)
        _cells[perm[i]].part = (i < half) ? 0 : 1;
}

void Partitioner::bfsPartition(mt19937& rng)
{
    int start = rng() % _cellNum;

    vector<bool> visited(_cellNum, false);
    vector<int> order;
    order.reserve(_cellNum);

    queue<int> q;
    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        int c = q.front(); q.pop();
        order.push_back(c);

        for (int n : _cellNets[c])
            for (int j : _netCells[n])
                if (!visited[j]) {
                    visited[j] = true;
                    q.push(j);
                }
    }

    for (int c = 0; c < _cellNum; ++c)
        if (!visited[c])
            order.push_back(c);

    int half = _cellNum / 2;
    _partSize[0] = half;
    _partSize[1] = _cellNum - half;
    for (int i = 0; i < _cellNum; ++i)
        _cells[order[i]].part = (i < half) ? 0 : 1;
}

void Partitioner::perturbPartition(mt19937& rng, int numSwaps)
{
    vector<int> g0, g1;
    g0.reserve(_partSize[0]);
    g1.reserve(_partSize[1]);
    for (int i = 0; i < _cellNum; ++i) {
        if (_cells[i].part == 0) g0.push_back(i);
        else                     g1.push_back(i);
    }
    shuffle(g0.begin(), g0.end(), rng);
    shuffle(g1.begin(), g1.end(), rng);

    int swaps = min(numSwaps, min((int)g0.size(), (int)g1.size()));
    for (int i = 0; i < swaps; ++i) {
        _cells[g0[i]].part = 1;
        _cells[g1[i]].part = 0;
    }
}

void Partitioner::saveBest()
{
    _bestCutSize = _cutSize;
    _bestPart.resize(_cellNum);
    for (int i = 0; i < _cellNum; ++i)
        _bestPart[i] = _cells[i].part;
}

void Partitioner::restoreBest()
{
    _cutSize = _bestCutSize;
    _partSize[0] = _partSize[1] = 0;
    for (int i = 0; i < _cellNum; ++i) {
        _cells[i].part = _bestPart[i];
        _partSize[_cells[i].part]++;
    }
}

void Partitioner::initNetDist()
{
    for (int n = 0; n < _netNum; ++n)
        _nets[n].dist[0] = _nets[n].dist[1] = 0;

    for (int c = 0; c < _cellNum; ++c) {
        int p = _cells[c].part;
        for (int n : _cellNets[c])
            _nets[n].dist[p]++;
    }
}

int Partitioner::computeCutSize() const
{
    int cut = 0;
    for (int n = 0; n < _netNum; ++n)
        if (_nets[n].dist[0] > 0 && _nets[n].dist[1] > 0)
            ++cut;
    return cut;
}

// ---------------------------------------------------------------------------
// [Opt 4] Merged gain computation + bucket building in one pass
// [Opt 2] Also initializes _unlockNum tracking
// ---------------------------------------------------------------------------

void Partitioner::initFmPass()
{
    int bs = bucketSize();
    for (int p = 0; p < 2; ++p) {
        _bHead[p].assign(bs, -1);
        _maxGIdx[p] = -1;
        _unlockNum[p] = 0;
    }

    for (int c = 0; c < _cellNum; ++c) {
        _cells[c].gain = 0;
        _cells[c].locked = false;
        _cells[c].bPrev = -1;
        _cells[c].bNext = -1;

        int F = _cells[c].part;
        int T = 1 - F;
        for (int n : _cellNets[c]) {
            if (_nets[n].dist[F] == 1) _cells[c].gain++;
            if (_nets[n].dist[T] == 0) _cells[c].gain--;
        }

        bucketInsert(F, c);
        _unlockNum[F]++;
    }
}

void Partitioner::bucketInsert(int part, int cid)
{
    int gi = gainIdx(_cells[cid].gain);
    int head = _bHead[part][gi];
    _bHead[part][gi] = cid;
    _cells[cid].bNext = head;
    _cells[cid].bPrev = -1;
    if (head != -1)
        _cells[head].bPrev = cid;
    if (gi > _maxGIdx[part])
        _maxGIdx[part] = gi;
}

void Partitioner::bucketRemove(int cid)
{
    int part = _cells[cid].part;
    int gi = gainIdx(_cells[cid].gain);
    int prev = _cells[cid].bPrev;
    int next = _cells[cid].bNext;

    if (prev == -1)
        _bHead[part][gi] = next;
    else
        _cells[prev].bNext = next;

    if (next != -1)
        _cells[next].bPrev = prev;

    if (_bHead[part][gi] == -1 && gi == _maxGIdx[part]) {
        while (_maxGIdx[part] >= 0 && _bHead[part][_maxGIdx[part]] == -1)
            --_maxGIdx[part];
    }
}

void Partitioner::updateGain(int cid, int delta)
{
    bucketRemove(cid);
    _cells[cid].gain += delta;
    bucketInsert(_cells[cid].part, cid);
}

// ---------------------------------------------------------------------------
// Cell move / undo
// ---------------------------------------------------------------------------

void Partitioner::doMove(int cid)
{
    int F = _cells[cid].part;
    int T = 1 - F;

    bucketRemove(cid);
    _cells[cid].locked = true;
    _unlockNum[F]--;

    for (int n : _cellNets[cid]) {
        int& dF = _nets[n].dist[F];
        int& dT = _nets[n].dist[T];

        // --- before distribution update ---
        if (dT == 0) {
            // net n is entirely in F; moving cid to T will cut it
            // all free neighbours in F see their TE shrink -> gain++
            for (int j : _netCells[n])
                if (j != cid && !_cells[j].locked)
                    updateGain(j, +1);
        } else if (dT == 1) {
            // the single cell in T will no longer be the sole representative
            for (int j : _netCells[n]) {
                if (_cells[j].part == T && !_cells[j].locked) {
                    updateGain(j, -1);
                    break;
                }
            }
        }

        --dF;
        ++dT;

        // --- after distribution update ---
        if (dF == 0) {
            // net n is now entirely in T -> all free cells on it see TE grow
            for (int j : _netCells[n])
                if (j != cid && !_cells[j].locked)
                    updateGain(j, -1);
        } else if (dF == 1) {
            // the remaining single cell in F now becomes sole representative
            for (int j : _netCells[n]) {
                if (j != cid && _cells[j].part == F && !_cells[j].locked) {
                    updateGain(j, +1);
                    break;
                }
            }
        }
    }

    _cells[cid].part = T;
    _partSize[F]--;
    _partSize[T]++;
}

void Partitioner::undoMove(int cid)
{
    int cur = _cells[cid].part;
    int orig = 1 - cur;
    _cells[cid].part = orig;
    _partSize[cur]--;
    _partSize[orig]++;
    for (int n : _cellNets[cid]) {
        _nets[n].dist[cur]--;
        _nets[n].dist[orig]++;
    }
}

// ---------------------------------------------------------------------------
// Single FM pass (returns true if cut size improved)
// ---------------------------------------------------------------------------

bool Partitioner::fmPass()
{
    initFmPass();

    _moveStack.clear();

    int accGain = 0;
    int maxAccGain = 0;
    int bestMoveCount = 0;
    int noImproveSteps = 0;
    int earlyStopThresh = max(50, _cellNum / 20);

    for (int step = 0; step < _cellNum; ++step) {
        bool can0 = _unlockNum[0] > 0 &&
                     (_partSize[0] - 1 >= _lowerBound) &&
                     (_partSize[1] + 1 <= _upperBound);
        bool can1 = _unlockNum[1] > 0 &&
                     (_partSize[1] - 1 >= _lowerBound) &&
                     (_partSize[0] + 1 <= _upperBound);

        int cand0 = (can0 && _maxGIdx[0] >= 0) ? _bHead[0][_maxGIdx[0]] : -1;
        int cand1 = (can1 && _maxGIdx[1] >= 0) ? _bHead[1][_maxGIdx[1]] : -1;

        if (cand0 == -1 && cand1 == -1) break;

        int toMove;
        if      (cand0 == -1) toMove = cand1;
        else if (cand1 == -1) toMove = cand0;
        else {
            int g0 = _cells[cand0].gain, g1 = _cells[cand1].gain;
            if      (g0 > g1) toMove = cand0;
            else if (g1 > g0) toMove = cand1;
            else {
                // [Opt 7] Tie-break: prefer smaller degree, then larger partition
                int d0 = (int)_cellNets[cand0].size();
                int d1 = (int)_cellNets[cand1].size();
                if (d0 != d1)
                    toMove = (d0 < d1) ? cand0 : cand1;
                else
                    toMove = (_partSize[0] >= _partSize[1]) ? cand0 : cand1;
            }
        }

        accGain += _cells[toMove].gain;
        _moveStack.push_back(toMove);
        doMove(toMove);

        if (accGain > maxAccGain) {
            maxAccGain = accGain;
            bestMoveCount = step + 1;
            noImproveSteps = 0;
        } else {
            ++noImproveSteps;
            if (noImproveSteps > earlyStopThresh) break;
            // [Opt 5] Impossible to recover if too far behind
            if (accGain < maxAccGain - _maxPinNum) break;
        }
    }

    for (int i = (int)_moveStack.size() - 1; i >= bestMoveCount; --i)
        undoMove(_moveStack[i]);

    if (maxAccGain > 0) {
        _cutSize -= maxAccGain;
        return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Top-level partition driver
// ---------------------------------------------------------------------------

void Partitioner::partition()
{
    _lowerBound = (int)ceil((1.0 - _bFactor) / 2.0 * _cellNum);
    _upperBound = (int)floor((1.0 + _bFactor) / 2.0 * _cellNum);

    auto tStart = chrono::steady_clock::now();
    auto elapsed = [&]() {
        return chrono::duration<double>(
            chrono::steady_clock::now() - tStart).count();
    };

    double timeBudget = min(280.0, max(5.0, _cellNum * 0.002));
    mt19937 rng(42);

    _moveStack.reserve(_cellNum);

    auto runFM = [&]() {
        initNetDist();
        _cutSize = computeCutSize();
        while (fmPass()) {
            if (elapsed() > timeBudget) break;
        }
    };

    // --- Phase 1: deterministic initial partition + FM ---
    initPartition();
    runFM();
    saveBest();

    // --- Phase 2: BFS-based initial partitions ---
    int numBFS = max(3, min(15, (int)(timeBudget * 0.15)));
    for (int i = 0; i < numBFS && elapsed() < timeBudget * 0.2; ++i) {
        bfsPartition(rng);
        runFM();
        if (_cutSize < _bestCutSize)
            saveBest();
    }

    // --- Phase 3: ILS - perturb best solution and refine ---
    int perturbSize = max(5, _cellNum / 500);
    int noImproveILS = 0;

    while (elapsed() < timeBudget * 0.90) {
        restoreBest();
        perturbPartition(rng, perturbSize);
        runFM();

        if (_cutSize < _bestCutSize) {
            saveBest();
            perturbSize = max(5, _cellNum / 500);
            noImproveILS = 0;
        } else {
            ++noImproveILS;
            if (noImproveILS % 3 == 0)
                perturbSize = min(perturbSize * 2, _cellNum / 4);
        }
    }

    // --- Phase 4: random restarts with remaining time ---
    while (elapsed() < timeBudget) {
        randomPartition(rng);
        runFM();
        if (_cutSize < _bestCutSize)
            saveBest();
    }

    restoreBest();
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
}

void Partitioner::writeResult(ofstream& outFile) const
{
    outFile << "Cutsize = " << _cutSize << '\n';

    outFile << "G1 " << _partSize[0] << '\n';
    for (int i = 0; i < _cellNum; ++i)
        if (_cells[i].part == 0)
            outFile << _cellName[i] << " ";
    outFile << ";\n";

    outFile << "G2 " << _partSize[1] << '\n';
    for (int i = 0; i < _cellNum; ++i)
        if (_cells[i].part == 1)
            outFile << _cellName[i] << " ";
    outFile << ";\n";
}
