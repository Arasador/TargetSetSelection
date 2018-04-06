#!/bin/bash

for file in create_instances/*; do
    #file = myciel2.col
    python maintopol.py 001000000 create_instances/${file##*/}
done
