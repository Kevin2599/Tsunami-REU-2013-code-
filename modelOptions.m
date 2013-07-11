function options = modelOptions(varargin)


	defaultOptions = struct('maxl',100,'timesteps',20000,'snapshots',100,'g',9.81, 'beach_slope',.05, 'bath_type', 'trapezoid', ...
                                 'dsigma',.01,'maxsigma',150,'xmax',5000,'DJN_beachwidth',50,'DJN_slopes',.5, ...
                            'trimAtBreak',false, ...
                            'timeFixStart',0.0, 'timeFixStride',10, 'timeFixEnd',0.1, ...
                            'quickLoad',false,'quickSave',false,'load',false,'save',false, ...
                            'plotLambda',true,'plotTime',false);
	savedOptions = struct();
	savedOptionsFile = '.REU_model_options';
    consoleOptions = readOptions(varargin);

  % load the savedOptions automatically (unless told otherwise)
	if ~readOption(consoleOptions,'dontLoadOptions',false)
		try
			load(savedOptionsFile);
		catch
		end
	end

  % save the userOptions so that they'll automatically load on subsequent runs
	userOptions = combineStructs(savedOptions,consoleOptions);
	if readOption(userOptions,'saveOptions',false)
    	userOptions.saveOptions = false;
    	savedOptions = userOptions;
    	save(savedOptionsFile,'savedOptions');
    end

	options = combineStructs(defaultOptions,userOptions);

    if ~isfield(options,'keeprate')
        options.keeprate = options.timesteps / options.snapshots;
    end

    options.bath.slope = options.beach_slope;
    if strcmp(options.bath_type,'trapezoid')
	    options.bath.trap_width = options.DJN_beachwidth;
	    options.bath.trap_slope = options.DJN_slopes;
	    options.bath.lr = @(h,bath) trapezoid_lr(h, bath.trap_width, bath.trap_slope);
	end

	if readOption(options,'printAllOptions',false)
		disp(options)
	end
end