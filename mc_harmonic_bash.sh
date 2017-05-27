#!/bin/bash

for k in {0..100};
do  
    ##g++ -std=c++17 mc_harmonic.cpp -o mc_harmonic
	i=$k  

	./mc_harmonic $i
done

./data_mc_harmonic

