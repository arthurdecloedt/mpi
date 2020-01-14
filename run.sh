#!/bin/zsh
rmtrash *.gol
mpicxx  -g -o  gol  -std=c++1z main.cpp
mpirun --oversubscribe --hostfile hostfile -np 4 gol 8 8 1 50 timings 1
mpirun --oversubscribe --hostfile hostfile -np 4 gol 10 10 1 50 timings
