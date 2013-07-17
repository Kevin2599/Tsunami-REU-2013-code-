%single value of lambda

function [phi psi] = evalPhi(func, lambda, sigma)
	% Phi = sum( (A.*sin(ln * lambda) + B.*cos(ln * l)) .* sinax_over_x(ln, sigma);

	A  = repmat( func.A ,1,length(sigma));
	B  = repmat( func.B ,1,length(sigma));
	[sigma ln] = meshgrid(sigma,func.lambda_n);

	phi = sum(ln .*  (A.*cos(ln .* lambda) - B.*sin(ln .* lambda)) .*   sinax_over_x(ln, sigma));
	psi = sum(       (A.*sin(ln .* lambda) + B.*cos(ln .* lambda)) .* d_sinax_over_x(ln, sigma));

	phi = phi(:); psi = psi(:);
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

	close0 = ax < .1;
	far0   = ~close0;

	y = zeros(size(ax));

	y(close0) = a(close0).^2 .* (-ax(close0)/3 + (ax(close0).^3)/30);
	y(far0  ) = (ax(far0) .* cos(ax(far0)) - sin(ax(far0))) ./ x(far0).^2;
end