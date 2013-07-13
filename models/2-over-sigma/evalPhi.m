function [phi psi] = evalPhi(func, lambda,sigma)
	% f_Phi = @(l,s) sum( (func.A.*sin(func.lambda_n * l) + func.B.*cos(func.lambda_n * l)) .* sin(func.lambda_n * s) / s);

	ln = func.lambda_n;
	A  = func.A;
	B  = func.B;


	f_phi = @(l,s) sum( ln .* (A.*cos(ln * l) - B.*sin(ln * l)) .*   sinax_over_x(ln, s));
	f_psi = @(l,s) sum(       (A.*sin(ln * l) + B.*cos(ln * l)) .* d_sinax_over_x(ln, s));


	phi = arrayfun(f_phi,lambda,sigma);
	psi = arrayfun(f_psi,lambda,sigma);
end

function y = d_sinax_over_x(a,x)
	% d/dx(f) ~ a^2 [-ax/3 + (ax)^3 / 30 - (ax)^5 / 840]
	% df < | (ax)^5 / 840 | / | -ax/3 + (ax)^3 / 30 |
	%    x<-2 (7 df+sqrt(7) sqrt(df (7 df+10)))
	%    with cutoff = .1, |f' - f|/f ~ 3.5e-7

	if(length(a) == 1)
		a = a * ones(size(x));
	elseif(length(x) == 1)
		x = x * ones(size(a));
	end

	ax = a.*x;

	close0 = ax <.1;
	far0   = ~close0;

	y = zeros(size(ax));

	y(close0) = a(close0).^2 .* (-ax(close0)/3 + (ax(close0).^3)/30);
	y(far0  ) = (ax(far0) .* cos(ax(far0)) - sin(ax(far0))) ./ x(far0).^2;
end