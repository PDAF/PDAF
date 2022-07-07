#!/bin/bash

cd offline_2D_serial
make clean
make cleandata
cd ..

cd offline_2D_parallel
make clean
make cleandata
cd ..

cd online_2D_serialmodel
make clean
make cleandata
cd ..

cd online_2D_serialmodel_2fields
make clean
make cleandata
cd ..

cd online_2D_parallelmodel
make clean
make cleandata
cd ..

cd online_2D_parallelmodel_fullpar
make clean
make cleandata
cd ..

cd online_2D_parallelmodel_fullpar_1fpe
make clean
make cleandata
cd ..

rm -f out.*
