function results = eval2Sig(options, bath)

	const = struct('sigma0', 2.9, 'dsigma', 0.01, 'maxsigma', 15, 'a', 0.1, 'p', 1, 'num_lambda', 40);
	lambda = 0:.1:10;

	println('Setting up initial conditions');

	sigma  = (0:options.dsigma:options.maxsigma)';
	lambda = (0:options.maxl/options.timesteps:options.maxl)';

	[phi_0 psi_0 ~] = exactPhi(sigma, zeros(size(sigma)), const);
	% phi_0 = phi_0 - mean(phi_0); % seems to make the initial condition fit better
	% psi_0 = zeros(size(sigma));

	println('Solving the equation');
	func = getPhi(sigma, phi_0, psi_0, options.num_lambda);

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

    results = struct('phi',phi', 'psi',psi', 'lambda',lambda ,'x0',DJN_x ,'eta0',DJN_eta, 'F',F, 'intF',intF);

	% figure(1); clf; hold on
	% plot(vars.t(2,:),vars.x(2,:),'b');
	% plot(exact.t(2,:),exact.x(2,:),'r');
	% figure(2); clf; hold on
	% x = abs(exact.x(2,:));
	% varsx = interp1(vars.t(2,:),vars.x(2,:),exact.t(2,:));
	% err = abs((exact.x(2,:) - varsx)./exact.x(2,:));
	% hold on
	% loglog(x,err,'b');
	% N = convhull(log(x),log(err));
	% loglog(x(N),err(N),'r');

end