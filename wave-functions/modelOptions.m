%{
    options = modelOptions(varargin)
		This loads options for the model to run.
		The passed options are read using 'readOptions'.
		These options are passed to plotWave, which has its own set of options.
		Option precedence: console (typed in by user), saved options, defaults.

    Model Options
	 -maxl [100]
	 -timesteps [20,000]
	 	lambda values that the model's evaluated at are 0:maxl/timesteps:maxl

	 -dsigma:		[.01]
	 -maxsigma:		[150]
		sigma values that the model's evaluated at are 0:dsigma:maxsigma

	 -g [9.81]
	 	gravity

	 -beach_slope [0.05]
	 	slope of the bay in the X direction (alpha)

	 -bath_type ['trapezoid']
	 	Defines the bathymetry of the bay on the yz plane

		>> trapezoid \__/ options
		 	-DJN_beachwidth [50]
		 	-DJN_slopes     [.5]
		 		The trapezoid's width and slope of the walls

	 -xmax [5000]
	 	?

	 -trimAtBreak [false]
		If the wave breaks, then the data after that point is removed

    Run Options
	 Data Retention
	  -snapshots [100]
	  -keeprate []
		The numerical model has ~100x more data points than necessary.
		These options trim the values of lambda for which data is saved.
		- Snapshots specifies how many values of lambda are saved.
		- Keeprate overrides snapshots and keeps every nth lambda value.

	  -timeFixStart  [.0]
	  -timeFixStride [10]
	  -timeFixEnd    [.1]
	  	Trims the values of sigma to (start*length:stride:end*length)

	 Save/Load
	  -save [false]
	  	Saves the data to external file for comparison with other models.

	  -quickSave [false]
	  -quickLoad [false]
	  	Caches the data from a run to make tweaking plots easier.

	 Plotting
	  -plotLambda [true]
	 	Plots x(lambda) eta(lambda) t(lambda)

	  -plotTime [false]
	 	Plots x(t) eta(t)

    Options Options
	 -saveOptions [false]
		saves the current saved/console options to file to be loaded automatically later

	 -dontLoadOptions
	 	prevents the saved options from being loaded
		combine with 'saveOptions' to clear the current saved options

	 -printAllOptions
		shows the value of all options currently defined

    See also: trapModel, readOptions
%}

function options = modelOptions(varargin)


	defaultOptions = struct('maxl',100,'timesteps',20000,'snapshots',100,'g',9.81, 'beach_slope',.05, 'model_type', 'trapezoid', ...
                                 'dsigma',.01,'maxsigma',150,'xmax',5000,'DJN_beachwidth',50,'DJN_slopes',.5, ...
                            'trimAtBreak',false, ...
                            'timeFixStart',0.0, 'timeFixStride',10, 'timeFixEnd',0.1, ...
                            'quickLoad',false,'quickSave',false,'save',false, ...
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
    if strcmp(options.model_type,'trapezoid')
    	options.model = @trapEval;
	    options.bath.trap_width = options.DJN_beachwidth;
	    options.bath.trap_slope = options.DJN_slopes;
	    options.bath.lr = @(h,bath) trapezoid_lr(h, bath.trap_width, bath.trap_slope);
	elseif strcmp(options.model_type,'2/sigma')

		options.model = @eval2Sig;
		options.bath.lr = @(h) [-(h.^.5) (h.^.5)];
		options.a = readOption(options,'a', .1);
		options.p = readOption(options,'p', 1);
		options.num_lambda = readOption(options,'num_lambda', 40);
		options.sigma0 = readOption(options,'sigma0', 2.9);
	end

	if readOption(options,'printAllOptions',false)
		disp(options)
	end
end