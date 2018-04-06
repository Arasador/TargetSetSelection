#!/bin/bash

for file in instances/*; do
    #file = myciel2.col
    python maintopol.py 1111111 instances/${file##*/} degree
done
