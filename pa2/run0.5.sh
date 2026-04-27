#!/bin/bash

input_folder="input_pa2"
evaluator_folder="evaluator"
output_folder="output"
output_file="summary_v1.txt"

# Overwrite the chosen summary file
: > "$output_file"

#run 
for input_name in ami33 ami49 apte hp xerox; do
    echo "---------REPORT FOR ${input_name}-----------" | tee -a "$output_file"

        echo "Running ./bin/fp 0.5 ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}.rpt" | tee -a "$output_file"
        ./bin/fp 0.5 ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}.rpt | tee -a "$output_file"
        ./${evaluator_folder}/evaluator.sh ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}.rpt 0.5 | tee -a "$output_file"
        echo "--------------------------------" | tee -a "$output_file"

    echo "--------------------------------" | tee -a "$output_file"
done