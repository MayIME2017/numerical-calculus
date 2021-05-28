function arg_curve = lagrangepoly(arg_nodes,arg_curve)
% Polynomial intrpolation with Lagrange method 
number_of_nodes = size(arg_nodes,1);
arg_curve(:,2)=0;
number_of_solutions=size(arg_curve,1);
pvals = zeros(number_of_nodes,number_of_nodes);
 
for k = 1 : number_of_solutions
    for j = 1 : number_of_nodes
        aux_product=1;
        for i = 1 : number_of_nodes
            if i ~= j
                aux_product = ...
                    aux_product*(arg_curve(k,1)-arg_nodes(i,1))/(arg_nodes(j,1)-arg_nodes(i,1));
                pp = poly(arg_nodes( (1:number_of_nodes) ~= i,1));
                pvals(i,:) = pp ./ polyval(pp, arg_nodes(i,1));
            end
        end
        arg_curve(k,2)=arg_curve(k,2)+aux_product*arg_nodes(j,2);
    end
end


end