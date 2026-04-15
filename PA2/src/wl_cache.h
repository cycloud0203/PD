#ifndef WL_CACHE_H
#define WL_CACHE_H

#include <climits>
#include "btree.h"
#include "shared_data.h"

class WLCache {
public:
    WLCache(int numBlocks, int numNets)
        : _numBlocks(numBlocks), _numNets(numNets) {}

    double compute(const BTree& tree, const SharedData& sd) const {
        const int* bx  = tree.blockXArr();
        const int* by  = tree.blockYArr();
        const int* bx2 = tree.blockX2Arr();
        const int* by2 = tree.blockY2Arr();

        long long totalWL2 = 0;
        for (int i = 0; i < _numNets; ++i) {
            int minX = INT_MAX, maxX = INT_MIN;
            int minY = INT_MAX, maxY = INT_MIN;
            for (int j = sd.netStart[i]; j < sd.netStart[i + 1]; ++j) {
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
            totalWL2 += (long long)(maxX - minX) + (maxY - minY);
        }
        return (double)totalWL2 * 0.5;
    }

    void accept() {}
    void invalidate() {}

private:
    int _numBlocks;
    int _numNets;
};

#endif
