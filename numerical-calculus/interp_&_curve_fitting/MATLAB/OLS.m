function [coefficients,arg_curve] = OLS(arg_nodes,arg_curve, polynomial_order)
% Calculates the polynomial least square fitting
number_of_nodes = size(arg_nodes,1);
A=zeros(number_of_nodes,polynomial_order+1);
% Constructing the matrix of powers of x of the data
for i = 1:number_of_nodes
    for j = 1:polynomial_order+1
        A(i,j)=arg_nodes(i,1)^(polynomial_order+1-j);
    end
end
% Solving the system Mu=y, where M is the matrix of powers of x and b
% u is the vector of the polynomial coefficients
% A'Aa=A'y - defining the square matrix A'A in order to solve the linear
% system
C=A'*A;
d=A'*arg_nodes(:,2);
coefficients=C\d;
arg_curve(:,2)=polyval(coefficients,arg_curve(:,1));
end