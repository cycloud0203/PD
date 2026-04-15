#!/bin/bash

input_folder="input_pa2"
evaluator_folder="evaluator"

for input_name in ami33 ami49 apte hp xerox; do
    echo "REPORT FOR ${input_name}" >> summary.txt
    for i in 0.1 0.3 0.5 0.7 0.9; do
        echo "Running ./bin/fp $i ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${input_folder}/${input_name}_${i}.rpt" >> summary.txt
        ./bin/fp $i ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${input_folder}/${input_name}_${i}.rpt >> summary.txt
        ./${evaluator_folder}/evaluator.sh ${input_folder}/${input_name}.block ${input_folder}/${input_name}.nets ${input_folder}/${input_name}_${i}.rpt $i >> summary.txt
        echo "--------------------------------" >> summary.txt
    done
    echo "--------------------------------" >> summary.txt
done