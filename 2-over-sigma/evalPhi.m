%single value of lambda

function [phi psi] = evalPhi(func, lambda, sigma)
	% Phi = sum( (A.*sin(ln * lambda) + B.*cos(ln * l)) .* sinax_over_x(ln, sigma);

	A  = repmat( func.A ,1,length(sigma));
	B  = repmat( func.B ,1,length(sigma));
	[sigma ln] = meshgrid(sigma,func.lambda_n);

	phi = sum(ln .*  (A.*cos(ln .* lambda) - B.*sin(ln .* lambda)) .* sinax_over_x(ln, sigma,0));
	psi = sum(       (A.*sin(ln .* lambda) + B.*cos(ln .* lambda)) .* sinax_over_x(ln, sigma,1));

	phi = phi(:); psi = psi(:);
end