set term pdfcairo enhanced color font "Arial,12" fontscale 0.5 size 12.7cm,8cm
set output 'compare.pdf'
set xlabel 'x'
set ylabel 'y'
set autoscale xy
set key top left

set style line 1 lt -1 lc rgb "blue" lw 2 dt 1
set style line 2 lt -1 lc rgb "purple" lw 2 dt 2
set style line 3 lt -1 lc rgb "red" lw 2 dt 3
set style line 4 lt -1 lc rgb "black" lw 2 dt 4
set style line 5 lt 2 lc rgb 'green' lw 1 dt 1
set style line 6 lt -1 lc rgb "yellow" lw 2 dt 10

set key ins vert
set key right top

set title "Comparison between lines of action"
plot 'trabalho4_FORTRAN_results.txt' using 1:2 with lines ls 1 title 'Fortran cubic spline',\
     'trabalho4_FORTRAN_results.txt' using 1:3 with lines ls 2 title 'Fortran approximation', \
     'trabalho4_FORTRAN_results.txt' using 1:4 with lines ls 3 title 'Fortran Lagrange interpolation',\
     'trabalho4_FORTRAN_results.txt' using 1:5 with lines ls 4 title 'Fortran Newton interpolation',\
     'trabalho4_nodes.txt' using 1:2 with points ls 5 title 'Nodes'

set title "Comparison between Cubic Spline codes"
plot 'trabalho4_FORTRAN_results.txt' using 1:2 with lines ls 1 title 'Fortran',\
     'trabalho4_FORTRAN_results.txt' using 1:6 with lines ls 2 title 'FGSL',\
     'trabalho4_spline_MATLAB_results.txt' using 1:2 with lines ls 3 title 'Matlab code',\
     'trabalho4_intrinsic_MATLAB_results.txt' using 1:2 with lines ls 4 title 'Matlab intrinsic',\
     'trabalho4_nodes.txt' using 1:2 with points ls 5 title 'Nodes'

set title "Comparison between Approximation codes"
plot 'trabalho4_FORTRAN_results.txt' using 1:3 with lines ls 1 title 'Fortran',\
     'trabalho4_approx_MATLAB_results.txt' using 1:2 with lines ls 3 title 'Matlab code',\
     'trabalho4_intrinsic_MATLAB_results.txt' using 1:3 with lines ls 4 title 'Matlab intrinsic',\
     'trabalho4_nodes.txt' using 1:2 with points ls 5 title 'Nodes',\

set title "Comparison between Polynomial Interpolation codes"
plot 'trabalho4_FORTRAN_results.txt' using 1:4 with lines ls 1 title 'Fortran',\
     'trabalho4_interpL_MATLAB_results.txt' using 1:2 with lines ls 3 title 'Matlab Lagrange code',\
     'trabalho4_interpN_MATLAB_results.txt' using 1:2 with lines ls 6 title 'Matlab Newton code',\
     'trabalho4_intrinsic_MATLAB_results.txt' using 1:4 with lines ls 4 title 'Matlab intrinsic (linear)',\
     'trabalho4_nodes.txt' using 1:2 with points ls 5 title 'Nodes',\
     'trabalho4_FORTRAN_results.txt' using 1:7 with lines ls 2 title 'FGSL',\

set term pop
