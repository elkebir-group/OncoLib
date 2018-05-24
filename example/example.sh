#!/bin/bash
build="../build"


echo "1. Simulating a clone tree corresponding to a metastatic tumor..."
$build/simulate -s 10 -m 1 -k 1 -kP 1 -C 160 -o .
echo "Output: T_seed10.tree T_seed10.vertex.labeling G_seed10.graph"
echo

echo "2. Visualizing simulated clone tree..."
$build/tree2dot -l T_seed10.vertex.labeling T_seed10.tree > T_seed10.dot
echo "Output: T_seed10.dot"
echo

echo "3. Visualizing simulated migration graph..."
$build/graph2dot G_seed10.graph > G_seed10.dot
echo "Output: G_seed10.dot"
echo

echo "4. Extracting mutation frequencies from clone tree..."
$build/tree2freqs T_seed10.tree > F_seed10.tsv
echo "Output: F_seed10.tsv"
echo

echo "5a. Generating samples as mixtures of clone tree leaves..."
$build/mix T_seed10.tree -k 2 > T_sampled_seed10.tree
echo "Output: T_sampled_seed10.tree"
echo
echo "5b. Visualizing resulting clone tree..."
$build/tree2dot T_sampled_seed10.tree > T_sampled_seed10.dot
echo "Output: T_sampled_seed10.dot"
echo

echo "6. Generating read matrix from sampled clone tree..."
$build/sequence -C 160 T_sampled_seed10.tree > R_sampled_seed10.tsv
echo "Output: R_sampled_seed10.tsv"
echo

echo "7. Generating frequency interval matrix from read matrix..."
$build/reads2freqs -a 0.05 R_sampled_seed10.tsv > F_sampled_seed10.tsv
echo "Output: F_sampled_seed10.tsv"
