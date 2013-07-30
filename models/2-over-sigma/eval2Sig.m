function results = eval2Sig(options, bath)
	println('Setting up initial conditions');

	sigma  = options.sigma;
	lambda = options.lambda;
	phi_0 = options.phi0;
	psi_0 = options.psi0;

	println('Solving the equation');
	func = getPhi(sigma, phi_0, psi_0, options.num_lambda);

	println('Evaluating at specific values')
	[lambda_mesh sigma_mesh] = meshgrid(lambda, sigma);
	figure(1); clf;
	% phi/psi = sigma x lambda
	for i=1:length(lambda)
		[phi1 psi1] = evalPhi(func, lambda(i), sigma);
		[e_phi e_psi e_Phi] = exactPhi(sigma, lambda(i)*ones(size(sigma)), const);
		 phi(:,i) =  phi1;  psi(:,i) =  psi1;
		ephi(:,i) = e_phi; epsi(:,i) = e_psi;
		plot(sigma, [phi1 e_phi phi1-e_phi]);
		title(sprintf('lambda (%.2f/%.2f)',lambda(i),lambda(end)));
		legend('phi','phi err');
		grid on
		drawnow();
	end


	println('Converting to physical variables')
	F    = (2/3) * sigma;
	intF = (1/3) * sigma.^2;

    results = struct('phi',phi', 'psi',psi', 'lambda',lambda, 'F',F, 'intF',intF);
end