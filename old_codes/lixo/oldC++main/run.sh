#!/bin/bash

make clean
make
./pci_test 0101011 preprocessed_instances/small_world-N1000-d4-p0.3-s1 >> out.txt
