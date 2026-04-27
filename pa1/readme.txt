Program: FM-based Hypergraph Partitioner

1. Build

Requirements:
- Linux environment
- g++ with C++17 support
- make

Compile command:
make

This builds the executable:
bin/fm

To remove the executable:
make clean

2. Run

Program usage:
./bin/fm <input_file> <output_file>

Example:
./bin/fm input_pa1/input_0.dat results/output_0.dat

3. Batch execution

To run all provided input cases (input_0.dat to input_5.dat):
bash run.sh

The script writes outputs to:
results/output_0.dat ... results/output_5.dat

4. Output format

The output file format is:
Cutsize = <value>
G1 <number_of_cells_in_G1>
<cell_name_1> <cell_name_2> ... ;
G2 <number_of_cells_in_G2>
<cell_name_1> <cell_name_2> ... ;

5. Evaluation

Checker and evaluator are provided in the evaluator directory.
Example checker usage:
./evaluator/checker input_pa1/input_0.dat results/output_0.dat

Example evaluator usage:
bash evaluator/evaluator.sh input_pa1/input_0.dat results/output_0.dat 0.10
