function modelOptions = sigmaOver2Options(userOptions)
	modelOptions = struct('num_lambda',40, 'model',@eval2Sig, 'maxsigma',15,'dlambda',.05, 'maxl',15, 'showError',false);
	modelOptions.phi0_psi0 = @(opt) sigmaOver2Init(opt);
	modelOptions.bath.lr = @(h) [-(h.^.5) (h.^.5)];
end