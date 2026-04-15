#ifndef SHARED_DATA_H
#define SHARED_DATA_H

#include <vector>

struct NetPin {
    int blockId;
    int fixedX2;
    int fixedY2;
};

struct SharedData {
    int numBlocks;
    int numNets;
    int outlineW, outlineH;
    double alpha;

    std::vector<int> blockW, blockH;
    std::vector<int> netStart;
    std::vector<NetPin> allPins;
    std::vector<std::vector<int>> blockNets;
    std::vector<std::vector<int>> blockAdj;

    double normA, normW;
};

#endif
