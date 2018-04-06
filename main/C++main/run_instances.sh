#!/bin/bash

for file in created_instances/*; do
	make clean
	make
	./pci_test 11111111 created_instances/${file##*/} >> out.txt
done