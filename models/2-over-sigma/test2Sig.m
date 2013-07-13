function [x_mesh t_samples eta_mesh] = test2Sig(const,lambda)
	println('Setting up initial conditions');

	sigma = (0:const.delSigma:const.maxSigma)';

	phi_0 = initialHump(const.a, const.sigma0, const.p, sigma);
	phi_0 = phi_0 - mean(phi_0);
	psi_0 = zeros(size(sigma));

	println('Solving the equation');
	func = getPhi(sigma, phi_0, psi_0, const.num_lambda);

	func.sigma = sigma;
	func.phi_0 = phi_0;

	println('Evaluating at specific values')
	[lambda_mesh sigma_mesh] = meshgrid(lambda, sigma);
	[phi psi] = evalPhi(func, lambda_mesh, sigma_mesh);
	figure(1); clf;
	for i=1:length(lambda)
		plot(sigma, phi(:,i));
		drawnow();
	end

	println('Converting to physical variables')
	F    = (2/3) * sigma;
	intF = (1/3) * sigma.^2; 
	options.g = 9.81; options.bath = struct('slope',.05);
	options.trimAtBreak = false;
	vars = convertToPhysicalVariables(phi', psi', lambda(:), F, intF, options);

	figure(1); clf; hold on
	for i=1:length(lambda)
		plot3(vars.x(:,i), vars.t(:,i), vars.eta(:,i));
		drawnow();
	end

	println('Reparameterizing from lambda to time')
	[x_mesh t_samples eta_mesh] = toConstantTime(vars.x(2:end-1,:), vars.t(2:end-1,:), 'linear', vars.eta(2:end-1,:));

	println('plotting')
	clf;
	for i=1:length(t_samples)
		plot(x_mesh(:,i), eta_mesh(:,i));
		xlabel('X'); ylabel('eta');
		title(['t = ' num2str(t_samples(i))]);
		drawnow();
	end
end

%{
	const = struct('sigma0', 2.9, 'delSigma', 0.01, 'maxSigma', 5, 'a', 0.1, 'p', 1, 'num_lambda', 400);
%}