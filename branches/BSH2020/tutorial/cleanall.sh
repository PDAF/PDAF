#!/bin/tcsh

cd offline_2D_serial
make clean
make cleandata
cd ..

cd offline_2D_openmp
make clean
make cleandata
cd ..

cd offline_2D_parallel
make clean
make cleandata
cd ..

cd online_2D_serialmodel_openmp
make clean
make cleandata
cd ..

cd online_2D_serialmodel_openmp_omi
make clean
make cleandata
cd ..

