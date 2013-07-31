function modelOptions = trapOptions(options)
	modelOptions = combineStructs( struct('DJN_beachwidth',50,'DJN_slopes',.5,'uFlag',0), options);

	modelOptions.model = @trapEval;
	modelOptions.bath.lr = @(h,bath) trapezoid_lr(h, bath.trap_width, bath.trap_slope);

	modelOptions.bath.trap_width = modelOptions.DJN_beachwidth;
	modelOptions.bath.trap_slope = modelOptions.DJN_slopes;

	modelOptions.initialWave = @trapInitialConditions;
end

function [x0 eta0 u0] = trapInitialConditions(options)
    %Define the initial profile
    x0 = -[0:1:options.xmax];

	%eta0=-9.0315e-4*exp(-1.5e-5*(1000+x).^2).*(1000+x);% N-wave \alpha =0.05
	% eta0=0.25 + 0.25*tanh((-x - 1000)/15);

 	eta0 = (-0.0001/0.6065) * exp(-2e-5 * (1000+x0).^2) .* (1000+x0); 

	eta0(abs(eta0)<1e-5) = 0;

    u0 = options.uFlag * eta0 .* sqrt(options.g./(-options.bath.slope*x0 - eta0));
end