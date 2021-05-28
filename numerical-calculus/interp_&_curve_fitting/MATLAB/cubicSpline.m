function [arg_curve,coefficients,polynomial_order] = cubicSpline(arg_nodes,arg_curve); 
 % Based on: Andrew Hastings (2021). cubicSpline (https://www.mathworks.com/
 % matlabcentral/fileexchange/5028-cubicspline), MATLAB Central File 
 % Exchange. Retrieved May 4, 2021.
 
number_of_breaks = size(arg_nodes,1);
number_of_gaps = number_of_breaks-1;

A = arg_nodes(:,2).';
B = zeros(1,number_of_breaks);
C = zeros(1,number_of_breaks);
D = zeros(1,number_of_breaks);
 
if number_of_breaks < 3
    % linear interpolation
    B(1) = (arg_nodes(2,2)-arg_nodes(1,2))/(arg_nodes(2,1)-arg_nodes(1,1));
    B(2) = B(1);
end

%
% step 1: preparation
%
D(1) = arg_nodes(2,1) - arg_nodes(1,1);
C(2) = (arg_nodes(2,2) - arg_nodes(1,2))/D(1);
for i = 2: number_of_gaps
  D(i) = arg_nodes(i+1,1) - arg_nodes(i,1);
  B(i) = 2.0*(D(i-1) + D(i));
  C(i+1) = (arg_nodes(i+1,2) - arg_nodes(i,2))/D(i);
  C(i) = C(i+1) - C(i);
end

%
% step 2: end conditions
%
B(1) = -D(1);
B(number_of_breaks) = -D(number_of_breaks-1);
C(1) = 0.0;
C(number_of_breaks) = 0.0;
if(number_of_breaks ~= 3)
  C(1) = C(3)/(arg_nodes(4,1)-arg_nodes(2,1)) - C(2)/(arg_nodes(3,1)-arg_nodes(1,1));
  C(number_of_breaks) = C(number_of_breaks-1)/(arg_nodes(number_of_breaks,1)-arg_nodes(number_of_breaks-2,1)) - C(number_of_breaks-2)/(arg_nodes(number_of_breaks-1,1)-arg_nodes(number_of_breaks-3,1));
  C(1) = C(1)*D(1)^2/(arg_nodes(4,1)-arg_nodes(1,1));
  C(number_of_breaks) = -C(number_of_breaks)*D(number_of_breaks-1)^2/(arg_nodes(number_of_breaks,1)-arg_nodes(number_of_breaks-3,1));
end

%
% step 3: forward elimination
%
for i = 2: number_of_breaks
  h = D(i-1)/B(i-1);
  B(i) = B(i) - h*D(i-1);
  C(i) = C(i) - h*C(i-1);
end

%
% step 4: back substitution
%
C(number_of_breaks) = C(number_of_breaks)/B(number_of_breaks);
for j = 1: number_of_gaps
  i = number_of_breaks-j;
  C(i) = (C(i) - D(i)*C(i+1))/B(i);
end

%
% step 5: compute spline coefficients
%
B(number_of_breaks) = (arg_nodes(number_of_breaks,2) - arg_nodes(number_of_gaps,2))/D(number_of_gaps) + D(number_of_gaps)*(C(number_of_gaps) + 2.0*C(number_of_breaks));
for i = 1: number_of_gaps
  B(i) = (arg_nodes(i+1,2) - arg_nodes(i,2))/D(i) - D(i)*(C(i+1) + 2.0*C(i));
  D(i) = (C(i+1) - C(i))/D(i);
  C(i) = 3.*C(i);
end
C(number_of_breaks) = 3.0*C(number_of_breaks);
D(number_of_breaks) = D(number_of_breaks-1);

coefficients(:,1) = A;
coefficients(:,2) = B;
coefficients(:,3) = C;
coefficients(:,4) = D;

polynomial_order = size(coefficients,2)-1;

% Use the above data A, B, C, and D, along with x to calculate the y spline points.
% if x is ouside the interval from spline definition take a boundary value
number_of_points = size(arg_curve,1);
for j = 1: number_of_gaps
    for i = 1:number_of_points
        if(arg_curve(i,1) <= arg_nodes(1,1))
            arg_curve(i,2) = arg_nodes(1,2);
        end
        if(arg_curve(i,1) >= arg_nodes(number_of_breaks,1))
            arg_curve(i,2) = arg_nodes(number_of_breaks,2);
        end
        %
        % search for arg_nodes(j,1) <= arg_curve(i,1) <= arg_nodes(j+1,1) 
        %
        if arg_curve(i,1) <= arg_nodes(j+1,1) && arg_curve(i,1) >= arg_nodes(j,1)
            dx = arg_curve(i,1) - arg_nodes(j,1);
            arg_curve(i,2) = A(j) + dx*(B(j) + dx*(C(j) + dx*D(j)));
        end
    end
end      
end