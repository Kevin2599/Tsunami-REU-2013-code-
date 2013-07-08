function [A B] = getPhi(sigma,g2eta0,U0,F,n)
	xeqtanx = @(ns) arrayfun((@(n) fzero((@(x) atan(x) + n*pi -x),[(n-.5)*pi (n+.5)*pi])),ns);
	lambda = xeqtanx(1:n);

	A_sol = g2eta0;
	B_sol = cumtrapz(sigma - [0 ;sigma(1:end-1)],U0 .* F);

	sin_lambda_sigma = sin(sigma * sqrt(lambda));
	A = ((1./sigma) * sqrt(lambda)) .* sin_lambda_sigma \ A_sol;
	B = ((1./sigma) * ones(size(lambda))) .* sin_lambda_sigma \ B_sol;
end