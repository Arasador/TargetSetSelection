#!/bin/bash
for file in coloring_problem_instances/*; do
	python functions_preprocess_instances/preprocess.py 10001 coloring_problem_instances/${file##*/} degree
done
