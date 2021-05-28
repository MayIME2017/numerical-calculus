function [coefficients, arg_curve] = poly_approx(arg_nodes,...
    arg_curve, polynomial_order)
% calculates the polynomial coefficients for each step
coefficients = polyfit(arg_nodes(:,1),arg_nodes(:,2),polynomial_order);
% Calculates the values for required points for each step
arg_curve(:,2) = polyval(coefficients, arg_curve(:,1));
end
