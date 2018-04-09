#!/bin/bash
make clean
make
for file in ../instances/small_world_instances/Instances5000/*
do
	./pci_test 1000000000000000 ../instances/small_world_instances/Instances5000/${file##*/}
done
