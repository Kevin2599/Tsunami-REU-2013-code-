function err = sigma2Err2(varargin)
	clf;


	err = [];
	sigma = (0:.01:15)';

	if false
		err_dsigma = .1;
		err_maxsigma = 10;

		sigma0s = 0:.3:15;
		p = 1;
		nlambdas = 1:40;
		for i = 1:length(sigma0s)
			sigma0 = sigma0s(i);

			[phi0 psi0] = sigmaOver2Init(struct('sigma0',sigma0,'a',1,'p',p,'sigma',sigma));
			int_phi0 = safeIntegrate(abs(phi0));
			for j = 1:length(nlambdas)
				nlambda = nlambdas(j);

				func = getPhi(sigma,phi0,psi0,nlambda);
				[phi1 psi1] = evalPhi(func,0,sigma);

				err(j,i) = safeIntegrate((phi1 - phi0).^2) / int_phi0;
				malprintf('\rsigma (%.2f/%.2f) n_lambda (%.2f/%.2f)',sigma0,sigma0s(end), nlambda, nlambdas(end));
			end
		end
		semilogy(err);
	end
	if true
		p = .05:.01:3;
		sigma0 = 0:.3:15;
		err = initDifference(struct('a',1,'nlambda',30,'sigma',sigma),'p',p,'sigma0',sigma0);

		[x y] = meshgrid(sigma0,p);
		mesh(x,y,err);
		xlabel('sigma0'); ylabel('p');
		set(gca, 'zscale', 'log')
		set(gca, 'zdir', 'reverse')
	end
end

function err = initDifference(defaults, var1_name,var1s, var2_name,var2s)
	for i=1:length(var1s)
		var1 = var1s(i);
		defaults = setfield(defaults,var1_name,var1);

		for j=1:length(var2s)
			var2 = var2s(j);
			defaults = setfield(defaults,var2_name,var2);


			[phi0 psi0] = sigmaOver2Init(defaults);
			func = getPhi(defaults.sigma, phi0,psi0, defaults.nlambda);
			[phi1 psi1] = evalPhi(func,0,defaults.sigma);

			err(i,j) = safeNorm(phi1 - phi0) / safeNorm(phi0);

			malprintf('\r%s (%.2f/%.2f) %s (%.2f/%.2f) total progress: %d/%d',var1_name, var1,var1s(end), var2_name,var2, var2s(end),...
				(i-1)*length(var2s) + j-1   ,    length(var2s)*length(var1s) );
		end
	end
	println('');
end

function iarr = safeNorm(arr, varargin)
	arr(isnan(arr)) = 0;
	iarr = norm(arr,varargin{:});
end