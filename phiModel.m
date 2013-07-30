 % phiModel(varagin)
	% This function takes a model which uses the Phi transform as described by Pelinovsky and Rybkin
	% It takes the results and:
	% 	Converts to physical variables
	% 	Aligns results with respect to t
	% 	Saves the results
	% 	Plots the results


function phiModel(varargin)
	options = modelSetup(varargin{:});
	bath = options.bath;

	if options.quickLoad
		load('.savedTrapModel.mat');
	else % run the model
		start_time = clock();

			println('Running model...')
		results = options.model(options);

			println('Converting to physical variables...')
		lambdaResults = convertToPhysicalVariables(results.phi, results.psi, results.lambda, results.F, results.intF, options);
		lambdaResults.lambda = results.lambda(:);


			println('Aligning with respect to time...')

		% get rid of lambda([1,end]) b/c they're NaN
		lambdaResults = applyFunToStruct(@(mat) mat(2:end-1,:) , lambdaResults);
		% doing the entire matrix is very slow
		lambdaResults = applyFunToStruct(@(mat) mat(floor(options.timeFixStart*end)+1 : options.timeFixStride : floor(options.timeFixEnd*end) ,:), lambdaResults);

		% z(x(l), t(l)) ... => z(x,t) ...
		[x_lin t_lin eta_lin u_lin] = alignToTime(lambdaResults.x, lambdaResults.t, 1:max(max(lambdaResults.t)) , lambdaResults.eta, lambdaResults.u);
		timeResults = struct('x',x_lin, 't',t_lin, 'eta',eta_lin, 'u',u_lin);

		fprintf('  - Simulation compeleted in %.2f seconds\n', etime(clock(),start_time));

	% Save the results
		if options.quickSave
			save('.savedTrapModel.mat',  'lambdaResults','timeResults','bath');
		end
		if options.save
			for i=1:size(x2,2)
				results.snapshot{i}.x=lambdaResults.x(:,i);
				results.snapshot{i}.eta=lambdaResults.eta(:,i);
				results.snapshot{i}.time=lambdaResults.t(2,i);
			end
			results.max_runup=max(max(lambdaResults.eta));
			results.case=['case_',num2str(bath.trap_width),'m_',num2str(1/bath.trap_slope),'_',num2str(bath.slope)];
			save(results.case,'results');
		end
	end

	println('Plotting...')
	if options.plotLambda
		plotWave(lambdaResults, bath, options);
		% plot(x0, eta0,'-b');
	end

	if options.plotTime
		plotWave( timeResults, bath, options, 'plotBathymetry',true, 'limitPlotT',100:140);
	end