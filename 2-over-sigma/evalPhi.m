function [phi phi_lam] = evalPhi(A,B,lambda_n, lambda,sigma)
	f_phi = @(l,s) sum( (A.*sin(lambda_n * l) + B.*cos(lambda_n * l)) .* sin(lambda_n * s) / s);

	d_lambda = .0001;

	phi = arrayfun(f_phi,lambda,sigma);
	phi_lam = (arrayfun(f_phi,lambda+d_lambda,sigma)-phi)./d_lambda;
end