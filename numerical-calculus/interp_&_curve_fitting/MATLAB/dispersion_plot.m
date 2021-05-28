function disp_fig_handle = dispersion_plot(arg_nodes)
disp_fig_handle = figure;
figure(disp_fig_handle)
plot(arg_nodes(:,1),arg_nodes(:,2),'kx');
xlabel('x')
ylabel('y')
end