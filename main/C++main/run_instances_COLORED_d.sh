#!/bin/bash

for file in ../instances/coloring_preprocessed_degree/*; do
	make clean
	make
	./pci_test 11111111 ../instances/coloring_preprocessed_degree/${file##*/}
done
