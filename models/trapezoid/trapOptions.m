function modelOptions = trapOptions(options)
	modelOptions = combineStructs( struct('DJN_beachwidth',50,'DJN_slopes',.5), options);

	modelOptions.model = @trapEval;
	modelOptions.bath.lr = @(h,bath) trapezoid_lr(h, bath.trap_width, bath.trap_slope);

	modelOptions.bath.trap_width = modelOptions.DJN_beachwidth;
	modelOptions.bath.trap_slope = modelOptions.DJN_slopes;
end