function plot2sig(x_mesh, t_samples, eta_mesh)
	clf;
	y_axis = [min(min(eta_mesh)) max(max(eta_mesh))];
	for i=1:length(t_samples);
		plot(x_mesh(:,i), eta_mesh(:,i));
		xlabel('X'); ylabel('eta');
		title(['t = ' num2str(t_samples(i))]);
		ylim(y_axis);
		drawnow();
	end
end