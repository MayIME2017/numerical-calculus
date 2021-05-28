function disp_fig_handle = gen_graphs(arg_nodes,arg_curve,arg_filename)
disp_fig_handle = figure('Visible', 'off');
figure(disp_fig_handle)
plot(arg_curve(:,1),arg_curve(:,2),'b-',arg_nodes(:,1),arg_nodes(:,2),'kx');
legend('Curve','Nodes')
xlabel('x')
ylabel('y')
print(disp_fig_handle, '-dpng', arg_filename);
end