#!/bin/bash

input_folder="input_pa2"
evaluator_folder="evaluator"
output_folder="output"

# if summary.txt exists, overwrite it
: > summary_full.txt

#run 
for input_name in ami33 ami49 apte hp xerox; do
    for i in 0.1 0.3 0.5 0.7 0.9; do
        echo "---------REPORT FOR ${input_name} _alpha = ${i}-----------" | tee -a summary.txt
        echo "Running ./bin/fp $i ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}_${i}.rpt" | tee -a summary_full.txt
        ./bin/fp $i ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}_${i}.rpt | tee -a summary_full.txt
        ./${evaluator_folder}/evaluator.sh ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${output_folder}/${input_name}_${i}.rpt $i | tee -a summary_full.txt
        echo "--------------------------------" | tee -a summary_full.txt
    done
done
echo "--------------------------------" | tee -a summary_full.txt