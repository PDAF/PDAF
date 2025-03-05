#!/bin/bash


cd offline_2D_parallel
make clean
make cleandata
cd ..

cd online_2D_parallelmodel
make clean
make cleandata
cd ..

rm -f out.*
