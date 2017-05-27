#!/bin/bash

lambda=5;

for k in {0..50..1};
do  
    ##g++ -std=c++17 mc_dot.cpp -o mc_dot
	i=$k  

	./mc_dot $i $lambda
done

./data_mc

echo "all done"

