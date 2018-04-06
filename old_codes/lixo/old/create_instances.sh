#!/bin/bash
for file in created_instances/*; do
	python preprocessed_instances/preprocess.py 10001 created_instances/${file##*/}
done