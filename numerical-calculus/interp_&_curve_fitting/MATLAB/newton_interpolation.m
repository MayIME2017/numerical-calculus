function arg_curve = newton_interpolation(arg_nodes,arg_curve)
% Polynomial interpolation with Divided Differences method
number_of_nodes = size(arg_nodes,1);
aux = zeros(1,number_of_nodes);
aux(1) = arg_nodes(1,2);
DivDiffTable = zeros(number_of_nodes-1,number_of_nodes-1);
for i = 1 : number_of_nodes - 1
   DivDiffTable(i,1) = (arg_nodes(i+1,2) - ...
       arg_nodes(i,2))/(arg_nodes(i+1,1) - arg_nodes(i,1));
end
for j = 2 : number_of_nodes - 1
   for i = 1 : number_of_nodes - j
      DivDiffTable(i, j) = (DivDiffTable(i+1, j - 1) - ...
          DivDiffTable(i, j - 1))/(arg_nodes(i+j,1) - arg_nodes(i,1));
   end
end
for j = 2 : number_of_nodes
   aux(j) = DivDiffTable(1, j-1);
end
product_x_parcels = ones(1,number_of_nodes);
DivDiff_parcels = aux(1);

for i = 1:size(arg_curve,1)
    for j = 2 : number_of_nodes
        product_x_parcels(j)=(arg_curve(i,1) - arg_nodes(j-1,1)) .* product_x_parcels(j-1);
        DivDiff_parcels(j) = aux(j) .* product_x_parcels(j);
    end
    arg_curve(i,2) = sum(DivDiff_parcels);
end
end