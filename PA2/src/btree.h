#ifndef BTREE_H
#define BTREE_H

#include <vector>
#include <random>
#include <utility>

struct BTreeNode {
    int parent;
    int left;
    int right;
    int blockId;
    bool rotated;
};

struct PerturbConfig {
    int rotateW = 7;
    int swapW = 6;
    int deleteInsertW = 6;
    int rebuildW = 1;
    int swapChildrenW = 4;
    int netAwarePercent = 50;
    int total() const { return rotateW + swapW + deleteInsertW + rebuildW + swapChildrenW; }

    static PerturbConfig forSize(int n) {
        if (n <= 12) return {6, 5, 5, 1, 5, 70};
        if (n <= 15) return {7, 6, 6, 1, 0, 0};
        if (n <= 50) return {6, 4, 8, 2, 0, 0};
        return {5, 2, 11, 2, 0, 0};
    }
};

class BTree
{
public:
    BTree(int numBlocks, const int* blockW, const int* blockH, std::mt19937& rng);

    void init();
    void randomInit();
    void randomInsertionInit();
    void netOrderedInit(const std::vector<std::vector<int>>& adjList);
    void pack();

    enum MoveType { ROTATE, SWAP, DELETE_INSERT, REBUILD, SWAP_CHILDREN };

    MoveType perturb(const PerturbConfig& cfg = PerturbConfig{});
    void undoPerturb();

    int getWidth() const { return _width; }
    int getHeight() const { return _height; }
    long long getArea() const { return (long long)_width * _height; }

    int getBlockX1(int blockId) const { return _blockX[blockId]; }
    int getBlockY1(int blockId) const { return _blockY[blockId]; }
    int getBlockX2(int blockId) const { return _blockX2[blockId]; }
    int getBlockY2(int blockId) const { return _blockY2[blockId]; }
    int getNumBlocks() const { return _numBlocks; }

    const int* blockXArr() const { return _blockX.data(); }
    const int* blockYArr() const { return _blockY.data(); }
    const int* blockX2Arr() const { return _blockX2.data(); }
    const int* blockY2Arr() const { return _blockY2.data(); }

    const std::vector<BTreeNode>& getNodes() const { return _nodes; }
    int getRoot() const { return _root; }
    void loadState(const std::vector<BTreeNode>& nodes, int root);

    void swapBlockIds(int nodeA, int nodeB);
    void rotateNode(int nodeIdx);
    void setRotation(int nodeIdx, bool rotated) { _nodes[nodeIdx].rotated = rotated; }
    void setBlockAdj(const std::vector<std::vector<int>>* adj) { _blockAdj = adj; }

private:
    std::vector<BTreeNode> _nodes;
    int _root;
    int _numBlocks;

    std::vector<int> _blockX, _blockY, _blockX2, _blockY2;
    int _width, _height;

    std::vector<int> _blockW;
    std::vector<int> _blockH;
    std::vector<int> _blockToNode;
    std::vector<int> _savedBlockToNode;

    std::vector<std::pair<int,int>> _contour;
    std::vector<std::pair<int,int>> _contourBuf;

    const std::vector<std::vector<int>>* _blockAdj = nullptr;
    std::mt19937& _rng;

    MoveType _lastMove;
    int _undoIdx1, _undoIdx2;
    bool _undoRot1;
    int _undoBid1, _undoBid2;
    bool _undoRot2;
    std::vector<BTreeNode> _savedNodes;
    int _savedRoot;

    void buildCompleteBinaryTree(const std::vector<int>& order, bool randRot);
    void deleteNode(int nodeIdx);
    void insertNode(int nodeIdx);
    void insertNodeNear(int nodeIdx, int targetNodeIdx);
};

#endif // BTREE_H
