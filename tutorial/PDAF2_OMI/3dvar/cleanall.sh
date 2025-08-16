#!/bin/tcsh

cd offline_2D_serial
make clean
make cleandata
cd ..

cd online_2D_serialmodel
make clean
make cleandata
cd ..

cd online_2D_parallelmodel
make clean
make cleandata
cd ..

rm -f out.*
