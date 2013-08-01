% function trapModel(options)
%
%   This program will find the solution to the tsunami runup problem on a
%   trapezoidal beach with a constant slope in the x direction.
%   It requires that the programs trapF.m and fixit.m be present.
%
%   Options are read in using 'readOptions'
%   For all options, see 'modelOptions'
%
%   The model is comprised of the functions:
%    trapF, setupModel, runModel, convertToPhysicalVariables, toConstantTime, plotWave
%
%   See also: readOptions, modelOptions

% The following varables are used in this program:
% W        - Vector that is used to find A.
% A        - n by n matrix that is used to solve the wave equation.
% a        - The amplitude of our gauss pulse.
% alpha    - The slope of the beach.
% b        - Length n vector that holds the right side of our system.
% dW       - The derivative of W. used to find A.
% dlambda  - The step size in lambda.
% dlambda2 - dlambda^2
% F        - Vector of length n that contains information about our cross
%            sections.
% g        - Gravity.
% intF     - Length n vector that is the integral of F. Use in conversion.
% intgrid  - Matrix used in the conversion.
% keeprate - Used in picking what values of lambda we will keep. Kept delta
%            lambda is 1/keeprate.
% l        - Counter to keep our information.
% lambda   - Contians out kept lambda values.
% maxl     - The maximum value of lambda.
% n        - The length of sigma.
% Phi_next - Length n vector that holds the next time step for our solution
%            to the wave equation.
% Phi_curr - Length n vector that holds the curent time step for our solution
%            to the wave equation. Needed to shuffle data.
% Phi_prev  - Length n vector that holds the previous time step for our solution
%            to the wave equation. Needed to shuffle data.
% Phiout   - n by length(lambda) matrix that contains out approxamation for
%            Phi.
% Psi_curr - Length n vector that holds the current time step for our solution
%            to the wave equation. Needed to shuffle data.
% Psi_prev - Length n vector that holds the previous time step for our solution
%            to the wave equation. Needed to shuffle data.
% Psiout   - n by length(lambda) matrix that contains out approxamation for
%            Psi.
% p        - The varience of our gauss pulse.
% s0       - The mean of our gauss pulse.
% sigma    - Vector that contains out values for sigma.
% step     - Counter that keeps track of lambda when solving our system
% timesteps- number of time steps between \lambda=0, and \lambda=maxl

function results = trapEval(options)

    timesteps = options.timesteps;            % number of time steps between \lambda=0, and \lambda=maxl, %DJN 4/10/13
    keeprate  = options.keeprate;     % keep every \it{keeprate}-th step.
    dlambda   = options.maxl/options.timesteps;
    dsigma    = options.dsigma;
    alpha     = options.bath.slope;
    g         = options.g;

    dsigma2  =  dsigma^2;
    dlambda2 = dlambda^2;

    println('Building model...')

  % Generate the bathymetry determined variables
    [sigma,F,H,H0,intF,dF,W,dW] = trapF(options, options.bath);
    W(1)=1e100; %W(1) is infinity, but just make it huge.

    n = length(sigma);

  % Holds the results for output
    Psiout = zeros(ceil(timesteps/keeprate), n);
    Phiout = zeros(ceil(timesteps/keeprate), n);
    lambda = zeros(ceil(timesteps/keeprate), 1);

  % Get the initial conditions
    [x0 eta0 u0] = options.initialWave(options);
    plot(x0, eta0)

    [Phi_prev Psi_prev] = convertToPhiPsi(x0,eta0,u0, g,alpha,H,F,sigma);
    % Phi_prev= -4*a*sigma.^(-1) .* ( (sigma-s0)/p^2 .* exp(-((sigma-s0)/p).^2) + (sigma+s0)/p^2 .* exp(-((sigma+s0)/p).^2));
    % Phi_prev(1)=0;

  % The matrix A solves for the next psi values in the system
    A=sparse(n,n); A(1,1) = 1; A(n,n) = 1;
    for i=2:n-1
        A(i, i-1) =   -(    dlambda2/(dsigma2) - dlambda2/(2*dsigma)*W(i)                  );
        A(i, i)   = 1 -( -2*dlambda2/(dsigma2)                            + dlambda2*dW(i) );
        A(i, i+1) =   -(    dlambda2/(dsigma2) + dlambda2/(2*dsigma)*W(i)                  );
    end
    % :Linear Boundary:
    % A(end,end)=dsigma+dlambda;
    % A(end,end-1)=-dlambda;



    println('Running model...')

  %record initial conditions
    Phiout(1,:) = Phi_prev;
    Psiout(1,:) = Psi_prev;
    lambda(1)=0;
    l = 2;

  % next step
    PHI_sigma = secondOrderDifference(Phi_prev,dsigma);
    PSI_sigma = secondOrderDifference(Psi_prev,dsigma) + W.*Psi_prev;

    Psi_curr = Psi_prev + dlambda*PHI_sigma;
    Phi_curr = Phi_prev + dlambda*PSI_sigma;

  % check if second step needs to be recorded
    if(keeprate == 1)
        Phiout(2,:) = Phi_curr;           
        Psiout(2,:) = Psi_curr;
        lambda(2)   = dlambda;
        l = 3;
    end

  % starts 2 steps in because the 2 previous steps are needed to calculate the next step
    for step=2:timesteps
         
        b = 2*Psi_curr - Psi_prev;
        b(1)=0; b(n)=0;

        % :Linear Boundary:
        % b(n)=dsigma*Psi_curr(n);
      
      % Solve for Psi
        Psi_prev = Psi_curr;
        Psi_curr = A\b;
        PSI_sigma = secondOrderDifference(Psi_curr,dsigma) + W.*Psi_curr;

      % Solve for Phi
        Phi_next = 4/3*Phi_curr - 1/3*Phi_prev + 2/3*PSI_sigma*dlambda;
        Phi_prev = Phi_curr;
        Phi_curr = Phi_next;
        
      % Record data at specific intervals
        if(mod(step,keeprate)==0)
            Psiout(l,:) = Psi_curr;
            Phiout(l,:) = Phi_curr;
            lambda(l)   = step*dlambda;

            if options.plotModelProgress
                hold off;
                plot(sigma(1:n-2),Phi_curr(1:n-2),'.b');
                set(gca,'xdir','reverse');
                title(['Step ' num2str(step) ' (' num2str(100 *step /options.timesteps) '%)']);
                drawnow();
            else
                printf('\r Step %d (%.1f%%)',step,100*step/options.timesteps);
            end

            l=l+1;
        end
    end
    println('');

    results = struct('phi',Phiout, 'psi',Psiout, 'lambda',lambda ,'x0',x0 ,'eta0',eta0, 'F',F, 'intF',intF);
end

function dv_dx = secondOrderDifference(v, dx)
    n = length(v);
    dv = zeros(n,1);

    dv(    1) = -3*v(  1) + 4*v(  2) - v(    3); % Second order  forward difference
    dv(2:n-1) =    v(3:n)            - v(1:n-2); % Second order  central difference
    dv(    n) = -3*v(  n) + 4*v(n-1) - v(  n-2); % Second order backward difference

    dv_dx = dv .* (1/(2*dx));
end



