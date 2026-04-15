#ifndef WL_CACHE_H
#define WL_CACHE_H

#include <climits>
#include <vector>
#include <cstdlib>
#include "btree.h"
#include "shared_data.h"

class WLCache {
public:
    WLCache(int numBlocks, int numNets)
        : _numBlocks(numBlocks), _numNets(numNets), _valid(false),
          _refCX2(numBlocks, 0), _refCY2(numBlocks, 0),
          _netWL(numNets, 0), _totalWL2(0), _pendingTotalWL2(0),
          _netDirty(numNets, 0) {}

    double compute(const BTree& tree, const SharedData& sd) {
        const int* bx  = tree.blockXArr();
        const int* by  = tree.blockYArr();
        const int* bx2 = tree.blockX2Arr();
        const int* by2 = tree.blockY2Arr();

        if (!_valid) {
            _totalWL2 = 0;
            for (int bid = 0; bid < _numBlocks; ++bid) {
                _refCX2[bid] = bx[bid] + bx2[bid];
                _refCY2[bid] = by[bid] + by2[bid];
            }
            for (int i = 0; i < _numNets; ++i) {
                _netWL[i] = computeNetWL(i, bx, by, bx2, by2, sd);
                _totalWL2 += _netWL[i];
            }
            _valid = true;
            _pendingTotalWL2 = _totalWL2;
            _movedBlocks.clear();
            _affectedNets.clear();
            return (double)_totalWL2 * 0.5;
        }

        _movedBlocks.clear();
        _movedCX2.clear();
        _movedCY2.clear();
        for (int bid = 0; bid < _numBlocks; ++bid) {
            int cx2 = bx[bid] + bx2[bid];
            int cy2 = by[bid] + by2[bid];
            if (cx2 != _refCX2[bid] || cy2 != _refCY2[bid]) {
                _movedBlocks.push_back(bid);
                _movedCX2.push_back(cx2);
                _movedCY2.push_back(cy2);
            }
        }

        if (_movedBlocks.empty()) {
            _pendingTotalWL2 = _totalWL2;
            _affectedNets.clear();
            return (double)_totalWL2 * 0.5;
        }

        _affectedNets.clear();
        _affectedNewWL.clear();
        for (int bid : _movedBlocks) {
            for (int netId : sd.blockNets[bid]) {
                if (!_netDirty[netId]) {
                    _netDirty[netId] = 1;
                    _affectedNets.push_back(netId);
                }
            }
        }

        _pendingTotalWL2 = _totalWL2;
        for (int netId : _affectedNets) {
            _netDirty[netId] = 0;
            _pendingTotalWL2 -= _netWL[netId];
            long long newWL = computeNetWL(netId, bx, by, bx2, by2, sd);
            _pendingTotalWL2 += newWL;
            _affectedNewWL.push_back(newWL);
        }

        return (double)_pendingTotalWL2 * 0.5;
    }

    void accept() {
        if (!_valid) return;
        _totalWL2 = _pendingTotalWL2;
        for (size_t i = 0; i < _movedBlocks.size(); ++i) {
            _refCX2[_movedBlocks[i]] = _movedCX2[i];
            _refCY2[_movedBlocks[i]] = _movedCY2[i];
        }
        for (size_t i = 0; i < _affectedNets.size(); ++i) {
            _netWL[_affectedNets[i]] = _affectedNewWL[i];
        }
    }

    void invalidate() { _valid = false; }

private:
    int _numBlocks;
    int _numNets;
    bool _valid;

    std::vector<int> _refCX2, _refCY2;
    std::vector<long long> _netWL;
    long long _totalWL2;
    long long _pendingTotalWL2;

    std::vector<char> _netDirty;
    std::vector<int> _movedBlocks;
    std::vector<int> _movedCX2, _movedCY2;
    std::vector<int> _affectedNets;
    std::vector<long long> _affectedNewWL;

    long long computeNetWL(int netId, const int* bx, const int* by,
                           const int* bx2, const int* by2,
                           const SharedData& sd) const {
        int minX = INT_MAX, maxX = INT_MIN;
        int minY = INT_MAX, maxY = INT_MIN;
        for (int j = sd.netStart[netId]; j < sd.netStart[netId + 1]; ++j) {
            int cx2, cy2;
            int bid = sd.allPins[j].blockId;
            if (bid >= 0) {
                cx2 = bx[bid] + bx2[bid];
                cy2 = by[bid] + by2[bid];
            } else {
                cx2 = sd.allPins[j].fixedX2;
                cy2 = sd.allPins[j].fixedY2;
            }
            if (cx2 < minX) minX = cx2;
            if (cx2 > maxX) maxX = cx2;
            if (cy2 < minY) minY = cy2;
            if (cy2 > maxY) maxY = cy2;
        }
        return (long long)(maxX - minX) + (maxY - minY);
    }
};

#endif
