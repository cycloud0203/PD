#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include "floorplanner.h"

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "Usage: ./fp <alpha> <input.block> <input.net> <output>" << std::endl;
        return 1;
    }

    double alpha = std::stod(argv[1]);

    std::fstream blockInput(argv[2], std::ios::in);
    if (!blockInput) {
        std::cerr << "Cannot open block file: " << argv[2] << std::endl;
        return 1;
    }

    std::fstream netInput(argv[3], std::ios::in);
    if (!netInput) {
        std::cerr << "Cannot open net file: " << argv[3] << std::endl;
        return 1;
    }

    std::fstream output(argv[4], std::ios::out);
    if (!output) {
        std::cerr << "Cannot open output file: " << argv[4] << std::endl;
        return 1;
    }

    auto startTime = std::chrono::steady_clock::now();

    Floorplanner fp(blockInput, netInput, alpha);
    fp.floorplan();

    auto endTime = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration<double>(endTime - startTime).count();

    fp.writeResult(output, runtime);

    blockInput.close();
    netInput.close();
    output.close();

    return 0;
}
