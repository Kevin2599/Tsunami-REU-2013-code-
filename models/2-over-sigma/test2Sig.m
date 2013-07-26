function [x_mesh t_samples eta_mesh] = test2Sig(const,lambda)
	if ~exist('const')
		const = struct('sigma0', 2.9, 'dsigma', 0.01, 'maxsigma', 15, 'a', 0.1, 'p', 1, 'num_lambda', 40);
		lambda = 0:.1:10;
	end

	println('Setting up initial conditions');

	sigma = (0:const.dsigma:const.maxsigma)';

	[phi_0 psi_0 ~] = exactPhi(sigma, zeros(size(sigma)), const);
	% phi_0 = phi_0 - mean(phi_0); % seems to make the initial condition fit better
	% psi_0 = zeros(size(sigma));

	println('Solving the equation');
	func = getPhi(sigma, phi_0, psi_0, const.num_lambda);

	func.sigma = sigma;
	func.phi_0 = phi_0;

	println('Evaluating at specific values')
	[lambda_mesh sigma_mesh] = meshgrid(lambda, sigma);
	figure(1); clf;
	% phi/psi = sigma x lambda
	for i=1:length(lambda)
		[phi1 psi1] = evalPhi(func, lambda(i), sigma);
		[e_phi e_psi e_Phi] = exactPhi(sigma, lambda(i)*ones(size(sigma)), const);
		 phi(:,i) =  phi1;  psi(:,i) =  psi1;
		ephi(:,i) = e_phi; epsi(:,i) = e_psi;
		plot(sigma, [phi1 phi1-e_phi]);
		title(sprintf('lambda (%.2f/%.2f)',lambda(i),lambda(end)));
		legend('phi','phi err');
		grid on
		drawnow();
	end


	println('Converting to physical variables')
	F    = (2/3) * sigma;
	intF = (1/3) * sigma.^2; 
	options.g = 9.81; options.bath = struct('slope',.05);
	options.trimAtBreak = false;
	vars = convertToPhysicalVariables(phi', psi', lambda(:), F, intF, options);
	exact = convertToPhysicalVariables(ephi',epsi', lambda(:), F, intF, options);

	figure(1); clf; hold on
	plot(vars.t(2,:),vars.x(2,:),'b');
	plot(exact.t(2,:),exact.x(2,:),'r');
	figure(2); clf; hold on
	x = abs(exact.x(2,:));
	varsx = interp1(vars.t(2,:),vars.x(2,:),exact.t(2,:));
	err = abs((exact.x(2,:) - varsx)./exact.x(2,:));
	hold on
	loglog(x,err,'b');
	N = convhull(log(x),log(err));
	loglog(x(N),err(N),'r');

	grid on;
	title('Sin series expansion of 2/sigma');
	ylabel('Relative error');
	xlabel('X value at water''s edge')
	return;


	figure(1); clf; hold on
	for i=1:length(lambda)
		plot3(vars.x(:,i), vars.t(:,i), vars.eta(:,i));
		drawnow();
	end

	println('Reparameterizing from lambda to time')
	[x_mesh t_samples eta_mesh] = toConstantTime(vars.x(2:end-1,:), vars.t(2:end-1,:), 'linear', vars.eta(2:end-1,:));

	println('plotting')
	% plotWave(struct('x',x_mesh, 't',t_samples, 'eta',eta_mesh), struct('slope',.05, 'lr', @(h) h.^.5));
end

%{
	const = struct('sigma0', 2.9, 'dsigma', 0.01, 'maxsigma', 15, 'a', 0.1, 'p', 1, 'num_lambda', 40);
%}