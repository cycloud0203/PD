#include <iostream>
#include <fstream>
#include "partitioner.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc != 3) {
        cerr << "Usage: ./fm <input file> <output file>" << endl;
        return 1;
    }

    ifstream input(argv[1]);
    if (!input) {
        cerr << "Cannot open input file: " << argv[1] << endl;
        return 1;
    }

    Partitioner partitioner;
    partitioner.parseInput(input);
    input.close();

    partitioner.partition();
    partitioner.printSummary();

    ofstream output(argv[2]);
    if (!output) {
        cerr << "Cannot open output file: " << argv[2] << endl;
        return 1;
    }
    partitioner.writeResult(output);
    output.close();

    return 0;
}
