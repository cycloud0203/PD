#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <string>
#include <vector>
#include <random>
#include <fstream>

struct Cell {
    int gain;
    int part;       // 0 = G1, 1 = G2
    bool locked;
    int bPrev;      // bucket doubly-linked list prev (-1 = head sentinel)
    int bNext;      // bucket doubly-linked list next (-1 = tail sentinel)
};

struct Net {
    int dist[2];    // number of cells in partition 0 and partition 1
};

class Partitioner {
public:
    void parseInput(std::ifstream& inFile);
    void partition();
    void writeResult(std::ofstream& outFile) const;
    void printSummary() const;

private:
    double _bFactor = 0;
    int _cellNum = 0;
    int _netNum = 0;
    int _maxPinNum = 0;
    int _cutSize = 0;
    int _partSize[2] = {0, 0};
    int _lowerBound = 0;
    int _upperBound = 0;

    std::vector<Cell> _cells;
    std::vector<std::string> _cellName;
    std::vector<std::vector<int>> _cellNets;

    std::vector<Net> _nets;
    std::vector<std::vector<int>> _netCells;

    // Bucket list per partition: _bHead[p][gainIdx] -> head cell id, -1 if empty
    std::vector<int> _bHead[2];
    int _maxGIdx[2] = {-1, -1};
    int _unlockNum[2] = {0, 0};

    // Persistent move stack (reused across FM passes to avoid reallocation)
    std::vector<int> _moveStack;

    // Best solution across all FM runs
    std::vector<int> _bestPart;
    int _bestCutSize = 0;

    // Gain index helpers
    int gainIdx(int gain) const { return gain + _maxPinNum; }
    int bucketSize() const { return 2 * _maxPinNum + 1; }

    void initPartition();
    void randomPartition(std::mt19937& rng);
    void bfsPartition(std::mt19937& rng);
    void perturbPartition(std::mt19937& rng, int numSwaps);
    void saveBest();
    void restoreBest();
    void initNetDist();
    void initFmPass();
    int  computeCutSize() const;

    void bucketInsert(int part, int cid);
    void bucketRemove(int cid);
    void updateGain(int cid, int delta);

    bool fmPass();
    void doMove(int cid);
    void undoMove(int cid);
};

#endif // PARTITIONER_H
