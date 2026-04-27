==============================================================
  Fixed-Outline Floorplanner  --  Compile & Run Instructions
==============================================================

1. Prerequisites
   - g++ with C++17 support (tested with g++ 11+)
   - GNU Make
   - POSIX threads (pthread)
   - **I run on EDA UNION U12**

2. Compile
   From the project root directory, run:

       make

   This produces the executable at  bin/fp .

   To rebuild from scratch:

       make clean && make

3. Run
   Usage:

       ./bin/fp <alpha> <input.block> <input.nets> <output.rpt>

   Arguments:
       alpha        Weight between area and wirelength (0.0 to 1.0).
       input.block  Path to the block description file.
       input.nets   Path to the net-list file.
       output.rpt   Path to the output report file.

   Example:

       ./bin/fp 0.5 input_pa2/ami49.block input_pa2/ami49.nets output/ami49.rpt

4. Batch Run
   A convenience script runs all five benchmarks with alpha = 0.5:

       bash run0.5.sh

   To run all benchmarks across multiple alpha values (0.1 to 0.9):

       bash run.sh

5. Evaluate
   After obtaining a report, use the provided evaluator:

       ./evaluator/evaluator.sh <input.block> <input.nets> <output.rpt> <alpha>

6. Visualize (optional, requires Python 3 + matplotlib)

       python3 scripts/visualize.py output/ami49.rpt \
           --block input_pa2/ami49.block \
           --nets  input_pa2/ami49.nets \
           --out   output/viz/ami49.png

   Or batch-visualize all reports under output/:

       python3 scripts/visualize.py --batch output/ --block-dir input_pa2/

7. Environment Variable
   FP_SEED  --  Set a fixed RNG seed for reproducible runs.
       Example:  FP_SEED=42 ./bin/fp 0.5 ...

8. Output Format (*.rpt)
   Line 1: cost  (alpha * area + (1 - alpha) * wirelength)
   Line 2: total wirelength
   Line 3: chip area
   Line 4: chip_width  chip_height
   Line 5: runtime (seconds)
   Line 6+: <block_name> <x1> <y1> <x2> <y2>
