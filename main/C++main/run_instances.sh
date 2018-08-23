#!/bin/bash
make clean
make

# Coloring problem instances
#for file in ../instances/coloring_problem_instances/*; do
#	./pci_test 00000001 ../instances/coloring_problem_instances/${file##*/} degree
#done

#./pci_test 00000001 ../instances/small_world_instances/Instances10000/small_world-N1000-d4-p0.3-s1
#./pci_test 00000001 ../instances/small_world_instances/Instances2500/small_world-N2500-d16-p0.3-s1
#./pci_test 00000001 ../instances/small_world_instances/Instances10000/small_world-N100-d4-p0.3-s1 small_world-N10000-d4-p0.3-s1

./pci_test 00000001 ../instances/small_world_instances/Instances10000/small_world-N10000-d4-p0.3-s1
