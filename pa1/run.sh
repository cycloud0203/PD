# ./fm input_pa1/input_0.dat results/output_0.dat
# bash evaluator.sh input_pa1/input_0.dat results/output_0.dat

# i from 0 to 5

for i in {0..5}; do
    ./bin/fm input_pa1/input_$i.dat results/output_$i.dat
done