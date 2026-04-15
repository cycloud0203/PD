#include "btree.h"
#include <algorithm>
#include <numeric>
#include <climits>

static const int NIL = -1;

BTree::BTree(int numBlocks, const int* blockW, const int* blockH, std::mt19937& rng)
    : _root(NIL), _numBlocks(numBlocks),
      _width(0), _height(0), _rng(rng),
      _lastMove(ROTATE), _undoIdx1(-1), _undoIdx2(-1),
      _undoRot1(false), _undoBid1(0), _undoBid2(0), _undoRot2(false),
      _savedRoot(NIL)
{
    _nodes.resize(_numBlocks);
    _blockX.resize(_numBlocks, 0);
    _blockY.resize(_numBlocks, 0);
    _blockX2.resize(_numBlocks, 0);
    _blockY2.resize(_numBlocks, 0);
    _contour.reserve(_numBlocks * 2 + 2);
    _contourBuf.reserve(_numBlocks * 2 + 2);
    _savedNodes.resize(_numBlocks);

    _blockW.assign(blockW, blockW + _numBlocks);
    _blockH.assign(blockH, blockH + _numBlocks);
    _blockToNode.resize(_numBlocks);
    _savedBlockToNode.resize(_numBlocks);

    init();
}

void BTree::buildCompleteBinaryTree(const std::vector<int>& order, bool randRot)
{
    for (int i = 0; i < _numBlocks; ++i) {
        _nodes[i].parent = NIL;
        _nodes[i].left = NIL;
        _nodes[i].right = NIL;
        _nodes[i].blockId = order[i];
        _nodes[i].rotated = randRot ? ((_rng() & 1) != 0) : false;
        _blockToNode[order[i]] = i;
    }
    if (_numBlocks == 0) return;
    _root = 0;
    for (int i = 0; i < _numBlocks; ++i) {
        int lc = 2 * i + 1, rc = 2 * i + 2;
        if (lc < _numBlocks) { _nodes[i].left = lc; _nodes[lc].parent = i; }
        if (rc < _numBlocks) { _nodes[i].right = rc; _nodes[rc].parent = i; }
    }
}

void BTree::init()
{
    std::vector<int> order(_numBlocks);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return (long long)_blockW[a] * _blockH[a] > (long long)_blockW[b] * _blockH[b];
    });
    buildCompleteBinaryTree(order, false);
}

void BTree::randomInit()
{
    std::vector<int> order(_numBlocks);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), _rng);
    buildCompleteBinaryTree(order, true);
}

void BTree::randomInsertionInit()
{
    std::vector<int> order(_numBlocks);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), _rng);

    for (int i = 0; i < _numBlocks; ++i) {
        _nodes[i].parent = NIL;
        _nodes[i].left = NIL;
        _nodes[i].right = NIL;
        _nodes[i].blockId = order[i];
        _nodes[i].rotated = (_rng() & 1) != 0;
        _blockToNode[order[i]] = i;
    }
    if (_numBlocks == 0) return;

    _root = 0;
    for (int i = 1; i < _numBlocks; ++i) {
        std::uniform_int_distribution<int> dist(0, i - 1);
        int target = dist(_rng);

        if (_rng() & 1) {
            int old = _nodes[target].left;
            _nodes[target].left = i;
            _nodes[i].parent = target;
            _nodes[i].left = old;
            _nodes[i].right = NIL;
            if (old != NIL) _nodes[old].parent = i;
        } else {
            int old = _nodes[target].right;
            _nodes[target].right = i;
            _nodes[i].parent = target;
            _nodes[i].left = NIL;
            _nodes[i].right = old;
            if (old != NIL) _nodes[old].parent = i;
        }
    }
}

void BTree::netOrderedInit(const std::vector<std::vector<int>>& adjList)
{
    if (_numBlocks == 0) return;

    int startBlock = 0;
    int maxDeg = -1;
    for (int i = 0; i < _numBlocks; ++i) {
        if ((int)adjList[i].size() > maxDeg) {
            maxDeg = (int)adjList[i].size();
            startBlock = i;
        }
    }

    std::vector<bool> visited(_numBlocks, false);
    std::vector<int> bfsOrder;
    std::vector<int> bfsParent;
    bfsOrder.reserve(_numBlocks);
    bfsParent.reserve(_numBlocks);

    bfsOrder.push_back(startBlock);
    bfsParent.push_back(-1);
    visited[startBlock] = true;
    int head = 0;
    while (head < (int)bfsOrder.size()) {
        int cur = bfsOrder[head];
        for (int nb : adjList[cur]) {
            if (!visited[nb]) {
                visited[nb] = true;
                bfsOrder.push_back(nb);
                bfsParent.push_back(head);
            }
        }
        head++;
    }
    for (int i = 0; i < _numBlocks; ++i) {
        if (!visited[i]) {
            bfsOrder.push_back(i);
            bfsParent.push_back(0);
        }
    }

    for (int i = 0; i < _numBlocks; ++i) {
        _nodes[i].parent = NIL;
        _nodes[i].left = NIL;
        _nodes[i].right = NIL;
        _nodes[i].blockId = bfsOrder[i];
        _nodes[i].rotated = (_rng() & 1) != 0;
        _blockToNode[bfsOrder[i]] = i;
    }

    _root = 0;
    for (int i = 1; i < _numBlocks; ++i) {
        int target;
        if (_rng() & 1)
            target = bfsParent[i];
        else {
            std::uniform_int_distribution<int> dist(0, i - 1);
            target = dist(_rng);
        }

        if (_rng() & 1) {
            int old = _nodes[target].left;
            _nodes[target].left = i;
            _nodes[i].parent = target;
            _nodes[i].left = old;
            _nodes[i].right = NIL;
            if (old != NIL) _nodes[old].parent = i;
        } else {
            int old = _nodes[target].right;
            _nodes[target].right = i;
            _nodes[i].parent = target;
            _nodes[i].left = NIL;
            _nodes[i].right = old;
            if (old != NIL) _nodes[old].parent = i;
        }
    }
}

void BTree::pack()
{
    _width = 0;
    _height = 0;
    if (_root == NIL) return;

    _contour.clear();
    _contour.push_back({0, 0});

    struct StackEntry { int nodeIdx; int parentX2; };
    StackEntry stkBuf[128];
    int stkTop = 0;
    stkBuf[stkTop++] = {_root, 0};

    while (stkTop > 0) {
        auto [nodeIdx, parentX2] = stkBuf[--stkTop];

        int bid = _nodes[nodeIdx].blockId;
        bool rot = _nodes[nodeIdx].rotated;
        int w = rot ? _blockH[bid] : _blockW[bid];
        int h = rot ? _blockW[bid] : _blockH[bid];

        int x1 = parentX2;
        int x2 = x1 + w;

        int maxY = 0;
        int yAtX2 = 0;
        int csz = (int)_contour.size();
        for (int i = 0; i < csz; ++i) {
            int sx = _contour[i].first;
            int sy = _contour[i].second;
            int nx = (i + 1 < csz) ? _contour[i + 1].first : INT_MAX;
            if (sx < x2 && nx > x1) {
                if (sy > maxY) maxY = sy;
            }
            if (sx <= x2 && nx > x2) {
                yAtX2 = sy;
            }
        }

        int y1 = maxY;
        int y2 = y1 + h;

        _blockX[bid] = x1;
        _blockY[bid] = y1;
        _blockX2[bid] = x2;
        _blockY2[bid] = y2;
        if (x2 > _width) _width = x2;
        if (y2 > _height) _height = y2;

        _contourBuf.clear();
        bool placed = false;
        bool x2done = (yAtX2 == y2);

        for (int i = 0; i < csz; ++i) {
            int sx = _contour[i].first;
            if (sx >= x1 && sx < x2) {
                if (!placed) {
                    _contourBuf.push_back({x1, y2});
                    placed = true;
                }
                continue;
            }
            if (!placed && sx >= x2) {
                _contourBuf.push_back({x1, y2});
                placed = true;
            }
            if (!x2done && sx >= x2) {
                if (sx == x2) x2done = true;
                else {
                    _contourBuf.push_back({x2, yAtX2});
                    x2done = true;
                }
            }
            _contourBuf.push_back(_contour[i]);
        }
        if (!placed) {
            _contourBuf.push_back({x1, y2});
        }
        if (!x2done) {
            _contourBuf.push_back({x2, yAtX2});
        }

        std::swap(_contour, _contourBuf);

        if (_nodes[nodeIdx].right != NIL) {
            stkBuf[stkTop++] = {_nodes[nodeIdx].right, x1};
        }
        if (_nodes[nodeIdx].left != NIL) {
            stkBuf[stkTop++] = {_nodes[nodeIdx].left, x2};
        }
    }
}

BTree::MoveType BTree::perturb(const PerturbConfig& cfg)
{
    int total = cfg.total();
    std::uniform_int_distribution<int> distOp(0, total - 1);
    int op = distOp(_rng);
    std::uniform_int_distribution<int> distNode(0, _numBlocks - 1);

    int rotEnd = cfg.rotateW;
    int swapEnd = rotEnd + cfg.swapW;
    int diEnd = swapEnd + cfg.deleteInsertW;
    int rebEnd = diEnd + cfg.rebuildW;

    if (op < rotEnd) {
        int idx = distNode(_rng);
        _undoIdx1 = idx;
        _nodes[idx].rotated = !_nodes[idx].rotated;
        _lastMove = ROTATE;
        return ROTATE;
    }
    else if (op < swapEnd) {
        int a = distNode(_rng);
        int b = distNode(_rng);
        while (b == a) b = distNode(_rng);
        _undoIdx1 = a;
        _undoIdx2 = b;
        _undoBid1 = _nodes[a].blockId;
        _undoRot1 = _nodes[a].rotated;
        _undoBid2 = _nodes[b].blockId;
        _undoRot2 = _nodes[b].rotated;
        std::swap(_nodes[a].blockId, _nodes[b].blockId);
        std::swap(_nodes[a].rotated, _nodes[b].rotated);
        _blockToNode[_nodes[a].blockId] = a;
        _blockToNode[_nodes[b].blockId] = b;
        _lastMove = SWAP;
        return SWAP;
    }
    else if (op < diEnd) {
        std::copy(_nodes.begin(), _nodes.end(), _savedNodes.begin());
        std::copy(_blockToNode.begin(), _blockToNode.end(), _savedBlockToNode.begin());
        _savedRoot = _root;
        int idx = distNode(_rng);
        int bid = _nodes[idx].blockId;
        deleteNode(idx);

        bool useNetAware = (cfg.netAwarePercent > 0
                            && _blockAdj != nullptr
                            && !(*_blockAdj)[bid].empty()
                            && ((int)(_rng() % 100) < cfg.netAwarePercent));
        if (useNetAware) {
            const auto& adj = (*_blockAdj)[bid];
            int k = std::min((int)adj.size(), 3);
            std::uniform_int_distribution<int> distK(0, k - 1);
            int neighborBid = adj[distK(_rng)];
            int targetNode = _blockToNode[neighborBid];
            insertNodeNear(idx, targetNode);
        } else {
            insertNode(idx);
        }
        _lastMove = DELETE_INSERT;
        return DELETE_INSERT;
    }
    else if (op < rebEnd) {
        std::copy(_nodes.begin(), _nodes.end(), _savedNodes.begin());
        std::copy(_blockToNode.begin(), _blockToNode.end(), _savedBlockToNode.begin());
        _savedRoot = _root;
        randomInit();
        _lastMove = REBUILD;
        return REBUILD;
    }
    else {
        int idx = distNode(_rng);
        _undoIdx1 = idx;
        std::swap(_nodes[idx].left, _nodes[idx].right);
        _lastMove = SWAP_CHILDREN;
        return SWAP_CHILDREN;
    }
}

void BTree::undoPerturb()
{
    switch (_lastMove) {
    case ROTATE:
        _nodes[_undoIdx1].rotated = !_nodes[_undoIdx1].rotated;
        break;
    case SWAP:
        _nodes[_undoIdx1].blockId = _undoBid1;
        _nodes[_undoIdx1].rotated = _undoRot1;
        _nodes[_undoIdx2].blockId = _undoBid2;
        _nodes[_undoIdx2].rotated = _undoRot2;
        _blockToNode[_undoBid1] = _undoIdx1;
        _blockToNode[_undoBid2] = _undoIdx2;
        break;
    case DELETE_INSERT:
    case REBUILD:
        std::copy(_savedNodes.begin(), _savedNodes.end(), _nodes.begin());
        std::copy(_savedBlockToNode.begin(), _savedBlockToNode.end(), _blockToNode.begin());
        _root = _savedRoot;
        break;
    case SWAP_CHILDREN:
        std::swap(_nodes[_undoIdx1].left, _nodes[_undoIdx1].right);
        break;
    }
}

void BTree::loadState(const std::vector<BTreeNode>& nodes, int root)
{
    std::copy(nodes.begin(), nodes.end(), _nodes.begin());
    _root = root;
    for (int i = 0; i < _numBlocks; ++i)
        _blockToNode[_nodes[i].blockId] = i;
}

void BTree::swapBlockIds(int a, int b)
{
    std::swap(_nodes[a].blockId, _nodes[b].blockId);
    std::swap(_nodes[a].rotated, _nodes[b].rotated);
}

void BTree::rotateNode(int idx)
{
    _nodes[idx].rotated = !_nodes[idx].rotated;
}

void BTree::deleteNode(int nodeIdx)
{
    int lc = _nodes[nodeIdx].left;
    int rc = _nodes[nodeIdx].right;
    int par = _nodes[nodeIdx].parent;

    int replacement = NIL;
    if (lc == NIL && rc == NIL) {
        replacement = NIL;
    } else if (lc == NIL) {
        replacement = rc;
    } else if (rc == NIL) {
        replacement = lc;
    } else {
        replacement = lc;
        int rightmost = lc;
        while (_nodes[rightmost].right != NIL)
            rightmost = _nodes[rightmost].right;
        _nodes[rightmost].right = rc;
        _nodes[rc].parent = rightmost;
    }

    if (replacement != NIL)
        _nodes[replacement].parent = par;

    if (par == NIL)
        _root = replacement;
    else if (_nodes[par].left == nodeIdx)
        _nodes[par].left = replacement;
    else
        _nodes[par].right = replacement;

    _nodes[nodeIdx].parent = NIL;
    _nodes[nodeIdx].left = NIL;
    _nodes[nodeIdx].right = NIL;
}

void BTree::insertNode(int nodeIdx)
{
    if (_root == NIL) {
        _root = nodeIdx;
        _nodes[nodeIdx].parent = NIL;
        _nodes[nodeIdx].left = NIL;
        _nodes[nodeIdx].right = NIL;
        return;
    }

    std::uniform_int_distribution<int> distNode(0, _numBlocks - 1);
    int target = distNode(_rng);
    while (target == nodeIdx) target = distNode(_rng);

    std::uniform_int_distribution<int> distSide(0, 1);
    if (distSide(_rng) == 0) {
        int oldChild = _nodes[target].left;
        _nodes[target].left = nodeIdx;
        _nodes[nodeIdx].parent = target;
        _nodes[nodeIdx].left = oldChild;
        _nodes[nodeIdx].right = NIL;
        if (oldChild != NIL) _nodes[oldChild].parent = nodeIdx;
    } else {
        int oldChild = _nodes[target].right;
        _nodes[target].right = nodeIdx;
        _nodes[nodeIdx].parent = target;
        _nodes[nodeIdx].left = NIL;
        _nodes[nodeIdx].right = oldChild;
        if (oldChild != NIL) _nodes[oldChild].parent = nodeIdx;
    }
}

void BTree::insertNodeNear(int nodeIdx, int targetNodeIdx)
{
    if (_root == NIL) {
        _root = nodeIdx;
        _nodes[nodeIdx].parent = NIL;
        _nodes[nodeIdx].left = NIL;
        _nodes[nodeIdx].right = NIL;
        return;
    }

    if (_rng() & 1) {
        int old = _nodes[targetNodeIdx].left;
        _nodes[targetNodeIdx].left = nodeIdx;
        _nodes[nodeIdx].parent = targetNodeIdx;
        _nodes[nodeIdx].left = old;
        _nodes[nodeIdx].right = NIL;
        if (old != NIL) _nodes[old].parent = nodeIdx;
    } else {
        int old = _nodes[targetNodeIdx].right;
        _nodes[targetNodeIdx].right = nodeIdx;
        _nodes[nodeIdx].parent = targetNodeIdx;
        _nodes[nodeIdx].left = NIL;
        _nodes[nodeIdx].right = old;
        if (old != NIL) _nodes[old].parent = nodeIdx;
    }
}
