#ifndef GLOBAL_BEST_H
#define GLOBAL_BEST_H

#include <mutex>
#include <vector>
#include "btree.h"

struct GlobalBest {
    std::mutex mtx;
    double bestCost = 1e18;
    bool feasible = false;
    std::vector<BTreeNode> nodes;
    int root = -1;
    int width = 0, height = 0;
    std::vector<int> bestX, bestY, bestX2, bestY2;

    bool tryUpdate(double cost, bool feas,
                   const std::vector<BTreeNode>& newNodes, int newRoot,
                   int w, int h,
                   const int* bx, const int* by,
                   const int* bx2, const int* by2,
                   int n)
    {
        std::lock_guard<std::mutex> lock(mtx);
        bool doUpdate = false;
        if (width == 0)                              doUpdate = true;
        else if (feas && !feasible)                  doUpdate = true;
        else if (feas == feasible && cost < bestCost) doUpdate = true;

        if (doUpdate) {
            bestCost = cost;
            feasible = feas;
            nodes    = newNodes;
            root     = newRoot;
            width    = w;
            height   = h;
            bestX.assign(bx, bx + n);
            bestY.assign(by, by + n);
            bestX2.assign(bx2, bx2 + n);
            bestY2.assign(by2, by2 + n);
        }
        return doUpdate;
    }
};

#endif
