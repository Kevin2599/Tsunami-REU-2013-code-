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

function results = trapEval(options)

    timesteps = options.timesteps;            % number of time steps between \lambda=0, and \lambda=maxl, %DJN 4/10/13
    keeprate  = options.keeprate;     % keep every \it{keeprate}-th step.
    dlambda   = options.maxl/options.timesteps;
    dsigma    = options.dsigma;
    alpha     = options.bath.slope;
    g         = options.g;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %We generate the space-determined variables sigma, F, H, H0, intF, dF, W,
    %and dW.

    [sigma,F,H,H0,intF,dF,W,dW] = trapF(options, options.bath);
    W(1)=1e100; %W(1) is infinity, but just make it huge.

    n = length(sigma);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build starting need information from user inputs and build the matrix A
    % that will be used to solve our system

    disp('Building model...')

    dsigma2  =  dsigma* dsigma;
    dlambda2 = dlambda*dlambda;

    A=sparse(n,n);
    b=zeros(n,1);


    A(1,1)=1;               % Define the matrix A, W and dW needed for our model
    for i=2:n-1
        A(i, i-1)=   -(    dlambda2/(dsigma2) - dlambda2/(2*dsigma)*W(i)                  );
        A(i, i)  = 1 -( -2*dlambda2/(dsigma2)                            + dlambda2*dW(i) );
        A(i, i+1)=   -(    dlambda2/(dsigma2) + dlambda2/(2*dsigma)*W(i)                  );
    end
    A(n,n)=1;


    %DJN
    %Define the initial profile
    DJN_x=-[0:1:options.xmax];
    %Comparison with the FUNWAVE model

    %eta_0 gives eta and u
    [DJN_eta,U_0]=eta_0(DJN_x);
    DJN_eta(abs(DJN_eta)<1e-5)=0;

    DJN_u=U_0*DJN_eta.*sqrt(g./(-alpha*DJN_x - DJN_eta));

    %DJN_eta=0.01*(1-tanh((1000+DJN_x)/200 ))/2
    plot(DJN_x, DJN_eta)

    %We need to convert (x, t, \eta, u) to (\sigma, \lambda, \phi, \psi)
    DJN_H=DJN_eta-DJN_x*alpha;
    DJN_Sigma=interp1(H, sigma, DJN_H);

    DJN_Phi=2*g*DJN_eta;


    %------------------MODEL STARTS HERE-----------------------
    %----------------------------------------------------------
    %----------------------------------------------------------
    % Define the initial Phi (wave height)
    % Phi_nm1=-4*a*sigma.^(-1).*((sigma-s0)/p^2.*exp(-1*((sigma-s0)/p).^2)+(sigma+s0)/p^2.*exp(-((sigma+s0)/p).^2));
    % Phi_nm1(1)=0;

    Phi_nm1=interp1(DJN_Sigma, DJN_Phi, sigma);
    Phi_nm1(isnan(Phi_nm1))=0;
    Phi_nm1(1)=0;
    Phi_nm1(end)=Phi_nm1(end-1);

    Phi_nm1=Phi_nm1(:);   %Make it a column, DJN 4/10/13

    PHI_sigma = secondOrderDifference(Phi_nm1)./(2*dsigma);

    DJN_u(isnan(DJN_u))=0;
    u_sigma=interp1(DJN_Sigma, DJN_u, sigma);
    u_sigma(isnan(u_sigma))=0;


    Psi_nm1=F.*u_sigma;

    Psi_nm1=Psi_nm1(:);   %Make it the column, DJN 4/10/13


    Psi_n=Psi_nm1+PHI_sigma*dlambda;                                     % Compute psi at the second step
    %DJN 4/10/13 %Psi=Psi_n;                                                   % Define Psi as the nth step



    PSI_sigma = secondOrderDifference(Psi_nm1)./(2*dsigma) + Psi_nm1.*W;

    Phi_n=Phi_nm1+dlambda*(PSI_sigma);                                                   % Compute phi at the nth step
    %DJN 4/10/13 %Phi=Phi_n;                                                                   % Define Psi as the nth step



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%FOR MOVING BOUNDRY OPTION%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the model for Psi and Phi

    disp('Running model...')

    %Pre-allocate for speed %DJN 4/10/13
    Psiout=zeros(ceil(timesteps/keeprate), n);
    Phiout=zeros(ceil(timesteps/keeprate), n);
    lambda=zeros(ceil(timesteps/keeprate), 1);

    step=0;
    l=1;                            % Index to keep only parts of our informaion
    Phiout(l,:)=Phi_nm1;         % Keep the initial conditions
    Psiout(l,:)=Psi_nm1;
    lambda(l)=0;

    %DJN 4/10/13, we need to keep the second step too.
    step=1;
    if(mod(step,keeprate)==0) %Check if we need to keep it.
        l=2;                            % Index to keep only parts of our informaion
        Phiout(l,:)=Phi_n;           
        Psiout(l,:)=Psi_n;
        lambda(l)=step*dlambda;
    end

    l=l+1;
    for step=2:timesteps    %we start from the third step, since the first two are already computed, DJN 4/10/13
    %DJN  b(1)=0;                     % Define b as the right side of our system
    %     for i=2:n-1
    %         b(i)=2*Psi_n(i)-Psi_nm1(i);
    %     end
    %     b(n)=0;
         
        b=2*Psi_n-Psi_nm1;          %Convert into the vector operation, DJN 4/10/13  
        b(1)=0; b(n)=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%Linear Boundary%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % A(end,end)=dsigma+dlambda;
        % A(end,end-1)=-dlambda;
        % b(n)=dsigma*Psi_n(n);
        
        
                
                
                
        Psi_nm1=Psi_n;              %We don't really need Psi vector, it just got eliminated to save time, DJN 4/10/13
        Psi_n=A\b;
        
        PSI_sigma = secondOrderDifference(Psi_n)./(2*dsigma) + W.*Psi_n;
        

        Phi = 4/3*Phi_n - 1/3*Phi_nm1 + 2/3*PSI_sigma*dlambda;                             % Define the next Phi
        Phi_nm1 = Phi_n;
        Phi_n = Phi;
          
        if(mod(step,keeprate)==0)              % Keep information at some points
            Psiout(l,:)=Psi_n;              % save the values at the current time step (written into the *_n arrays)
            Phiout(l,:)=Phi_n;
            lambda(l)=step*dlambda;

            hold off;
            plot(sigma(1,1:n-2),Phiout(l,1:n-2),'.b')
            set(gca,'xdir','reverse')
            title(['Step ' num2str(step) ' (' num2str(100 *step /options.timesteps) '%)'])
            drawnow();

            l=l+1;
        end
    end

    results = struct('phi',Phiout, 'psi',Psiout, 'lambda',lambda ,'x0',DJN_x ,'eta0',DJN_eta, 'F',F, 'intF',intF);

    % save(['analytical_nw_',results.case,'tmp.mat'], 'results')
end

function d = secondOrderDifference(v)
    n = length(v);
    d = zeros(length(v),1);
    d(    1) = -3*v(  1) + 4*v(  2) - v(    3); % Second order  forward difference
    d(2:n-1) =    v(3:n)            - v(1:n-2); % Second order  central difference
    d(    n) = -3*v(  n) + 4*v(n-1) - v(  n-2); % Second order backward difference
end