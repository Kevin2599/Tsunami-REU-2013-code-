function [phi psi Phi] = exactPhi(sigma,lambda,vars)
	d = 1e-6;
	exact_Phi = @(sigma,lambda) exact__Phi(sigma,lambda,vars);

	Phi = exact_Phi(sigma, lambda);
	phi = (exact_Phi(sigma,lambda+d) - Phi)./d;
	psi = (exact_Phi(sigma+d,lambda) - Phi)./d;
end

function Phi = exact__Phi(sigma,lambda,vars)
	e = @(lamSign, sig0Sign) gaussianFunction(sig0Sign.*vars.sigma0, vars.p, sigma + lamSign*lambda);
	Phi = (vars.a ./ sigma) .* ( e(1,-1) - e(-1,-1) + e(1,1) - e(-1,1) );
end