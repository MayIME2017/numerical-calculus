#!/bin/bash

make -f Makefile
rm -f *.o *.mod *~ *.exe
time ./trabalho4_fortran

gnuplot plot.plt

./trabalho4_value
