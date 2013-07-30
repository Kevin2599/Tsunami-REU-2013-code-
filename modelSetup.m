%{
    options = modelSetup(varargin)
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

function options = modelSetup(varargin)


	defaultOptions = struct('maxl',100,'timesteps',20000,'snapshots',100,'g',9.81,'dsigma',.01,'maxsigma',150,'xmax',5000, ...
							'beach_slope',.05, 'model_type', 'trapezoid', ...
							'a',.1, 'p',1, 'sigma0',2.9, ...
                            'trimAtBreak',false, ...
                            'timeFixStart',0.0, 'timeFixStride',10, 'timeFixEnd',0.1, ...
                            'quickLoad',false,'quickSave',true,'save',false, ...
                            'plotLambda',true,'plotTime',false);
	
	defaultOptions.phi0_psi0 = @(opt) deal(opt.a * gaussianFunction(opt.sigma0, opt.p, opt.sigma,1), zeros(size(opt.sigma)));

	savedOptions = struct();
	savedOptionsFile = '.REU_model_options';
    consoleOptions = readOptions(varargin);

  % load the savedOptions automatically (unless told otherwise)
	if ~readOption(consoleOptions,'dontLoadOptions',false)
		try; load(savedOptionsFile); catch; end
	end

  % save the userOptions so that they'll automatically load on subsequent runs
	userOptions = combineStructs(savedOptions,consoleOptions);
	if readOption(consoleOptions,'saveOptions',false)
    	userOptions.saveOptions = false;
    	savedOptions = userOptions;
    	save(savedOptionsFile,'savedOptions');
    end
  
  	modelOptions = struct();
  	userOptions.model_type = readOption(userOptions,'model_type','trapezoid');
    if strcmp(userOptions.model_type,'trapezoid')
    	modelOptions = trapOptions(userOptions);
	elseif strcmp(userOptions.model_type,'2/sigma')
		modelOptions = sigmaOver2Options(userOptions);
	else
		error(['Model ''' userOptions.model_type ''' not found']);
	end
	options = combineStructs(defaultOptions,modelOptions,userOptions);

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Setup from parameters

    if ~isfield(options,'keeprate')
        options.keeprate = options.timesteps / options.snapshots;
    end
    if ~isfield(options,'dlambda')
    	options.dlambda = options.maxl/options.timesteps;
    end

    options.bath.slope = options.beach_slope;
    options.lambda = (0:options.dlambda:options.maxl)';
	options.sigma  = (0:options.dsigma :options.maxsigma)';

  % Model Specific parameters / defaults  

    if isfield(options,'x0')
    	error('Conversion from initial x,eta,u is Not Implemented Yet.')
    elseif ~isfield(options,'phi0')
    	[phi_0 psi_0] = options.phi0_psi0(options);
    	options.phi0 = phi_0; options.psi0 = psi_0;
    end

	if readOption(options,'printAllOptions',false)
		disp(options)
	end
end