#!/bin/bash

for a in {1..10}
do
  python small_world_random_graph.py 2500 16 0.3 $a
done
#python small_world_random_graph.py <N>5000 <density>8 <p>0.3 <seed>i
