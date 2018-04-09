#!/bin/bash

for file in ../instances/coloring_problem_instances/*; do
	make clean
	make
	./pci_test 11111111 ../instances/coloring_problem_instances/${file##*/} degree
done
