% function trapModel(options)
%
%	This program will find the solution to the tsunami runup problem on a
%	trapezoidal beach with a constant slope in the x direction.
%	It requires that the programs trapF.m and fixit.m be present.
%
%	Options are read in using 'readOptions'
%	For all options, see 'modelOptions'
%
%	The model is comprised of the functions:
%	 trapF, setupModel, runModel, convertToPhysicalVariables, toConstantTime, plotWave
%
% SEE ALSO
% readOptions, modelOptions

% The following varables are used in this program:
% W        - Vector that is used to find A.
% A        - n by n matrix that is used to solve the wave equation.
% I        - Matrix to look for breaks in time.
% a        - The amplitude of our gauss pulse.
% alpha    - The slope of the beach.
% b        - Length n vector that holds the right side of our system.
% breakc   - Checks to see if we have broken at that time.
% brokeat  - Keeps the index if there was a break.
% dW       - The derivative of W. used to find A.
% dlambda  - The step size in lambda.
% dlambda2 - dlambda^2
% dummy    - Matrix that is not used but comes from building a grid.
% eta1     - n by length(lambda) matrix that contains our exact solution if
%            it exist.
% eta2     - n by length(lambda) matrix that contains our numerical
%            solution.
% Exact    - Bool that it true if a exact analytical solution is known to
%            exist.
% F        - Vector of length n that contains information about our cross
%            sections.
% Fgrid    - Matrix that is used to convert from nonphysical varables to
%            physical ones.
% G        - Length n vector that is used to solve the wave equation.
% g        - Gravity.
% i        - Counter.
% intF     - Length n vector that is the integral of F. Use in conversion.
% intgrid  - Matrix used in the conversion.
% keeprate - Used in picking what values of lambda we will keep. Kept delta
%            lambda is 1/keeprate.
% LAM      - Matrix that  is used to convert from nonphysical varables to
%            physical ones.
% l        - Counter to keep our information.
% lambda   - Contians out kept lambda values.
% leg      - Used to move legend.
% maxl     - The maximum value of lambda.
% n        - The length of sigma.
% Phi      - Length n vector that holds the curent time step for our solution
%            to the wave equation.
% Phi_n    - Length n vector that holds the curent time step for our solution
%            to the wave equation. Needed to shuffle data.
% Phi_nm1  - Length n vector that holds the curent time step for our solution
%            to the wave equation. Needed to shuffle data.
% Phiout   - n by length(lambda) matrix that contains out approxamation for
%            Phi.
% plotb    - bool to turn on plot.
% Psi      - Length n vector that holds the curent time step for our solution
%            to the wave equation.
% Psi_n    - Length n vector that holds the current time step for our solution
%            to the wave equation. Needed to shuffle data.
% Psi_nm1  - Length n vector that holds the current time step for our solution
%            to the wave equation. Needed to shuffle data.
% Psiout   - n by length(lambda) matrix that contains out approxamation for
%            Psi.
% p        - The varence of our gauss pulse.
% phi      - Exact phi if it exist.
% psi      - Exact psi if it exist.
% slope    - Finds the slope of the wave to check for breaking
% s0       - The mean of our gauss pulse.
% SIG      - Matrix that  is used to convert from nonphysical varables to
%            physical ones.
% sigma    - Vector that contains out values for sigma.
% step     - Counter that keeps track of lambda when solving our system
% t1       - Our time output for the exact solution. NOTE MATRIX
% t2       - Our time output for the aprox solution. NOTE MATRIX
% timesteps- number of time steps between \lambda=0, and \lambda=maxl %DJN 4/10/13
% u1       - n by length(lambda) matrix for our velocity output fot the
%            exact solution.
% u2       - n by length(lambda) matrix for our velocity output fot the
%            aprox solution.
% x1       - Our distance output for the exact solution. NOTE MATRIX
% x2       - Our distance output for the aprox solution. NOTE MATRIX

function trapModel(varargin)

	options = modelOptions(varargin{:});

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define all needed user inputs

	bath = options.bath;

	if options.quickLoad
		load('.savedTrapModel.mat');
	else % run the model
		start_time = clock();
		println('Building bathymetry...')

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% We generate the space-determined variables sigma, F, H, H0, intF, dF, W, and dW.
		[sigma,F,H,H0,intF,dF,W,dW] = trapF(options, 2*(bath.trap_width) + 2*bath.slope * (options.xmax/bath.trap_slope) );
		W(1)=1e100; %W(1) is the infinity, just make it huge, instead of the Inf, DJN 4/10/13
		W = W';
		
		%For no potential.
		%W=0*W;
		%dW=0*dW;
		
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Build the starting conditions from the user inputs
		println('Building model...')
		[Phi_nm1,Phi_n, Psi_nm1,Psi_n,counter,A,dlambda,W,PHI_LAMBDA, DJN_x, DJN_eta] = setupModel(sigma,W,dW,H,F,options);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Solve the model for Psi and Phi
		println('Running model...')
		[Phiout Psiout lambda] = runModel(sigma,Phi_nm1,Phi_n,Psi_nm1,Psi_n,counter,A,dlambda,W,PHI_LAMBDA,F,options);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Convert back to physical varables
		println('Converting Approx data...')
		lambdaResults = convertToPhysicalVariables(Phiout, Psiout, lambda, F, intF, options);
		lambdaResults.lambda = lambda;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Align data with respect to time (currently with respect to lambda)
		println('Aligning with respect to time...')

		% get rid of lambda=[1,end] b/c they fuck griddata up
		lambdaResults = applyFunToStruct(@(mat) mat(2:end-1,:) , lambdaResults);
		% doing the entire matrix is very slow
		lambdaResults = applyFunToStruct(@(mat) mat(floor(options.timeFixStart*end)+1 : options.timeFixStride : floor(options.timeFixEnd*end) ,:), lambdaResults);

		% z(x(l), t(l)) ... => z(x,t) ...
		[x_lin t_lin eta_lin u_lin] = toConstantTime(lambdaResults.x, lambdaResults.t, 1:max(max(lambdaResults.t)) , lambdaResults.eta, lambdaResults.u);
		timeResults = struct('x',x_lin, 't',t_lin, 'eta',eta_lin, 'u',u_lin);


	% Save the results
		fprintf('  - Simulation compeleted in %d seconds\n', ceil(etime(clock(),start_time)));

		if options.quickSave
			save('.savedTrapModel.mat',  'lambdaResults','timeResults','bath', 'DJN_x','DJN_eta');
		end
		if options.save
			for i=1:size(x2,2)
				results.snapshot{i}.x=x2(:,i);
				results.snapshot{i}.eta=eta2(:,i);
				results.snapshot{i}.time=t2(2,i);
			end
			results.max_runup=max(max(eta2));
			results.case=['case_',num2str(bath.trap_width),'m_',num2str(1/bath.trap_slope),'_',num2str(bath.slope)];
			save(results.case,'results');
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the data
%%%%%

	println('Plotting...')



	x_axis   = [min(min(timeResults.x))  , max(max(timeResults.x))+10];
	eta_axis = [min(min(timeResults.eta)), max(max(timeResults.eta)) ];

	max_height = eta_axis(2) - x_axis(1)*bath.slope;
	max_y      = max_height/bath.trap_slope + bath.trap_width/2;

	bath.height = [max_height; 0];
	bath.left   = [-max_y; -bath.trap_width/2];
	bath.right  = [ max_y;  bath.trap_width/2];

	% x=-3*max(max(x2)):.1:2*max(max(x2));
	if options.plotLambda
		plotWave(lambdaResults, bath, options);
		plot(DJN_x, DJN_eta,'-b');
	end

	if options.plotTime
		plotWave( timeResults, bath, options, 'plotBathymetry',true, 'limitPlotT',100:140);
	end

end %function
