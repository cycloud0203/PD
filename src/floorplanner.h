#ifndef FLOORPLANNER_H
#define FLOORPLANNER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <chrono>
#include "module.h"
#include "shared_data.h"
#include "global_best.h"

class Floorplanner
{
public:
    Floorplanner(std::fstream& blockInput, std::fstream& netInput, double alpha);
    ~Floorplanner();

    void floorplan();
    void writeResult(std::fstream& output, double runtime);

private:
    void parseBlock(std::fstream& input);
    void parseNet(std::fstream& input);
    void buildNetArrays();
    void buildSharedData();
    void computeNormalization();
    double computeWirelength();

    double _alpha;
    int _outlineW, _outlineH;

    std::vector<Block*> _blocks;
    std::vector<Terminal*> _terminals;
    std::vector<Net*> _nets;
    std::map<std::string, Terminal*> _nameToTerm;

    int _numNets;
    std::vector<int> _netStart;
    std::vector<NetPin> _allPins;
    std::vector<std::vector<int>> _blockNets;
    std::vector<std::vector<int>> _blockAdj;

    SharedData _sd;
    GlobalBest _gb;

    std::mt19937 _rng;
};

#endif // FLOORPLANNER_H
