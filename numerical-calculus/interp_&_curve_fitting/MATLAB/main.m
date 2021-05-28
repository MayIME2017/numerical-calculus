for i = 1:2
    points_spline=points;
    points_approx=points;
    points_interapprox = points;
    switch(i)
        case 1
            % Calculates the cubic spline interpolation
            [points_spline,coefficients_spline,polynomial_order_spline] = cubicSpline(nodes,points_spline);
        
            % Calculates the ordinary least squares approximation
            [coefficients_approx, points_approx] = ...
                OLS(nodes,points_approx,polynomial_order_approx);
        
            % Calculates the polynomial interpolation of approximation points
            nodes_interapprox(:,1) = nodes(:,1);
            nodes_interapprox(:,2) = polyval(coefficients_approx,nodes(:,1));
        
            % Calculates the Lagrange polynomial interpolation for the approximation
            % data
            points_interapprox = lagrangepoly(nodes_interapprox,points_interapprox);
            
            results(:,1) = points_spline(:,1);
            results(:,2) = points_spline(:,2);
            results(:,3) = points_approx(:,2);
            results(:,4) = points_interapprox(:,2);
        
            % Calculates the Divided Difference polynomial for the approximation
            % data
            points_interapprox = newton_interpolation(nodes_interapprox,points_interapprox);

            results(:,5) = points_interapprox(:,2);
        
        case 2
            % calculates the cubic spline interpolation
            pp = spline(nodes(:,1),nodes(:,2));
            [breaks,coefficients_spline,number_of_pieces,order,target] = ...
                unmkpp(pp);
            points_spline(:,2) = ppval(pp,points_spline(:,1));
            polynomial_order_spline = size(coefficients_spline,2)-1;
        
            % Calculates the ordinary least squares approximation for the data
            % The function defined uses polyfit and polyval from MATLAB
            % documentation - theses functions utilize the least-square method
            % to calculate the polynomial with the defined order. This means
            % that the chosen functions are of the form fk=akx^k
            [coefficients_approx, points_approx] = ...
                poly_approx(nodes,points_approx,polynomial_order_approx);
            % Calculates the interpolation polynomial from the approximation
            % points
            nodes_interapprox(:,1) = nodes(:,1);
            nodes_interapprox(:,2) = polyval(coefficients_approx,nodes(:,1));
            points_interapprox(:,2) = interp1(nodes_interapprox(:,1),nodes_interapprox(:,2),points_interapprox, 'linear');
            results(:,6) = points_spline(:,2);
            results(:,7) = points_approx(:,2);
            results(:,8) = points_interapprox(:,2);
    end
end

results_spline(:,1) = results(:,1);
results_spline(:,2) = results(:,2);
results_approx(:,1) = results(:,1);
results_approx(:,2) = results(:,3);
results_interpL(:,1) = results(:,1);
results_interpL(:,2) = results(:,4);
results_interpN(:,1) = results(:,1);
results_interpN(:,2) = results(:,5);
results_intrinsic(:,1) = results(:,1);
results_intrinsic(:,2) = results(:,6);
results_intrinsic(:,3) = results(:,7);
results_intrinsic(:,4) = results(:,8);