function func = getPhi(sigma,phi0,psi0,n)
	xeqtanx = @(ns) arrayfun((@(n) fzero((@(x) atan(x) + n*pi -x),[(n-.5)*pi (n+.5)*pi])),ns);
	lambda = xeqtanx(0:n) ./ sigma(end);

	if sigma(1) == 0
		phi0 = phi0(2:end);
		psi0 = psi0(2:end);
		sigma = sigma(2:end);
	end

	% B_sol = cumtrapz(sigma_step * ones(size(U0)), U0 .* F);

	sin_lambda_sigma = sin(sigma * lambda);

	func.A = ((1./sigma) * lambda) .* sin_lambda_sigma \ phi0(:);
	func.B = ((1./sigma) * ones(size(lambda))) .* sin_lambda_sigma \ psi0(:);

	func.lambda_n = lambda';

	% condition = cond(((1./sigma) * lambda) .* sin_lambda_sigma);
	% lambda = lambda(:);
	% difference = sum(abs(A_sol - evalPhi(A,B,lambda,0,sigma')));
end