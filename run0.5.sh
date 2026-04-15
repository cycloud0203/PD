#!/bin/bash

input_folder="input_pa2"
evaluator_folder="evaluator"
output_folder="output"

# if summary.txt exists, overwrite it
: > summary.txt

#run 
for input_name in ami33 ami49 apte hp xerox; do
    echo "---------REPORT FOR ${input_name}-----------" | tee -a summary.txt

        echo "Running ./bin/fp $i ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}_${i}.rpt" | tee -a summary.txt
        ./bin/fp 0.5 ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}_${i}.rpt | tee -a summary.txt
        ./${evaluator_folder}/evaluator.sh ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}_${i}.rpt 0.5 | tee -a summary.txt
        echo "--------------------------------" | tee -a summary.txt

    echo "--------------------------------" | tee -a summary.txt