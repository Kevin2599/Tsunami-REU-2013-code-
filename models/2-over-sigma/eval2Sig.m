function results = eval2Sig(options, bath)
	println('Setting up model...');

	sigma  = options.sigma;
	lambda = options.lambda;

	func = getPhi(sigma, options.phi0, options.psi0, options.num_lambda);

	println('Evaluating at specific values...')
	figure(1); clf;
	% phi/psi = sigma x lambda
	for i=1:length(lambda)
		[phi1 psi1] = evalPhi(func, lambda(i), sigma);
		 phi(:,i) =  phi1;  psi(:,i) =  psi1;

		if options.showError
			[e_phi e_psi e_Phi] = exactPhi(sigma, lambda(i)*ones(size(sigma)), const);
			ephi(:,i) = e_phi; epsi(:,i) = e_psi;

			plot(sigma, [phi1 phi1-e_phi]);
			legend('phi','phi err');
		else
			plot(sigma, phi1);
		end
		title(sprintf('lambda (%.2f/%.2f)',lambda(i),lambda(end)));
		grid on
		drawnow();
	end

	F    = (2/3) * sigma;
	intF = (1/3) * sigma.^2;

    results = struct('phi',phi', 'psi',psi', 'lambda',lambda, 'F',F, 'intF',intF);
end