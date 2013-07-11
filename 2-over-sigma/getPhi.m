function [A B lambda condition difference] = getPhi(L,eta0,U0,F,n)
	xeqtanx = @(ns) arrayfun((@(n) fzero((@(x) atan(x) + n*pi -x),[(n-.5)*pi (n+.5)*pi])),ns);
	lambda = xeqtanx(1:n)./L;
	%lambda = sqrt(lambda); % everytime after this, lambda_n is only seen as sqrt(lambda_n). This makes things simpler
	g = 9.81;

	eta0 = eta0(2:end);
	U0 = U0(2:end);
	F = F(2:end);

	sigma_step = L/length(F);
	sigma = (sigma_step:sigma_step:L)';

	A_sol = 2*g*eta0;
	B_sol = cumtrapz(sigma_step * ones(size(U0)), U0 .* F);

	sin_lambda_sigma = sin(sigma * lambda);
	A = ((1./sigma) * lambda) .* sin_lambda_sigma \ A_sol(:);
	B = ((1./sigma) * ones(size(lambda))) .* sin_lambda_sigma \ B_sol(:);

% hold off
% 	plot(A_sol,'r');
% hold on
% 	plot((((1./sigma) * lambda) .* sin_lambda_sigma) * A);

	condition = cond(((1./sigma) * lambda) .* sin_lambda_sigma);
	lambda = lambda(:);
	difference = sum(abs(A_sol - evalPhi(A,B,lambda,0,sigma')));
end