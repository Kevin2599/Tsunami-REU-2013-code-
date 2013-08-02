function [phi0 psi0] = sigmaOver2Init(vars)
	e = @(sigSign, sig0Sign) gaussianFunction(sig0Sign.*vars.sigma0, vars.p, sigSign * vars.sigma,1);
	phi0 = (vars.a ./ vars.sigma) .* ( e(1,-1) - e(-1,-1) + e(1,1) - e(-1,1) );
	psi0 = zeros(size(phi0));
end