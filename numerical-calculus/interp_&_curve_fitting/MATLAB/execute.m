clear all
close all
clc

problem_name = 'trabalho4';
problem_data_file = strcat(problem_name, '_data.m');

tic
run(problem_data_file)
run('main.m')
time = toc;

%SAVING RESULTS AND DATA POINTS TO PLAIN BINARY FILES
problem_nodes_file = strcat(problem_name,'_nodes.txt');
save(problem_nodes_file, 'nodes', '-ascii');
problem_results_file = strcat(problem_name,'_spline','_MATLAB_results.txt');
save(problem_results_file, 'results_spline', '-ascii');
problem_results_file = strcat(problem_name,'_approx','_MATLAB_results.txt');
save(problem_results_file, 'results_approx', '-ascii');
problem_results_file = strcat(problem_name,'_interpL','_MATLAB_results.txt');
save(problem_results_file, 'results_interpL', '-ascii');
problem_results_file = strcat(problem_name,'_interpN','_MATLAB_results.txt');
save(problem_results_file, 'results_interpN', '-ascii');
problem_results_file = strcat(problem_name,'_intrinsic','_MATLAB_results.txt');
save(problem_results_file, 'results_intrinsic', '-ascii');