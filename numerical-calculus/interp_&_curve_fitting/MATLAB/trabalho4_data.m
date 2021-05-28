%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT DATA FOR INTERPOLATION AND APPROXIMATION %
% ---------------------------------------------- %
% Written  by Mayara Magalhães Carvalho          %
% Instituto Militar de Engenharia - IME          %
% Last modified in May 2021                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes used to construct the function
nodes(:,1) = [1.59, 4.46, 7.32, 10.29, 13.31, 16.28, 19.33, 22.40, 25.29, ...
     28.17, 31.21, 34.32, 37.22, 40.20, 43.25, 46.25, 49.20];
nodes(:,2) = [92.75, 79.16, 76.56, 67.39, 65.92, 61.00, 50.76, 49.86, ...
    43.83, 40.74, 38.44, 34.87, 30.39, 28.70, 25.88, 22.72, 21.48];

% Number of nodes for the interpolation and approximation
number_of_nodes = size(nodes,1);

% Points to calculate the value from the built function
points(:,1) = 2:1:40;

% Polynomial order for the approximation
polynomial_order_approx = 3;