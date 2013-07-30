function modelOptions = sigmaOver2Options(userOptions)
	modelOptions = struct('num_lambda',40, 'model',@eval2Sig, 'maxsigma',15,'dlambda',.05, 'maxl',15);
	modelOptions.phi0_psi0 = @(opt) sigmaOver2Init(opt);
	modelOptions.bath.lr = @(h) [-(h.^.5) (h.^.5)];
end
function [phi0 psi0] = sigmaOver2Init(vars)
	e = @(sigSign, sig0Sign) gaussianFunction(sig0Sign.*vars.sigma0, vars.p, sigSign * vars.sigma,1);
	phi0 = (vars.a ./ vars.sigma) .* ( e(1,-1) - e(-1,-1) + e(1,1) - e(-1,1) );
	psi0 = zeros(size(phi0));
end