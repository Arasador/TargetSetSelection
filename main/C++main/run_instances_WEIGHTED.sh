#!/bin/bash

for file in created_instances/created_instances_weigthed/*; do
	make clean
	make
	./pci_test 11111111 created_instances/created_instances_weigthed/${file##*/}
done
