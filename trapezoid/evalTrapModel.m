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
%	See also: readOptions, modelOptions

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

function results = evalTrapModel(options)
	bath = options.bath;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% We generate the space-determined variables sigma, F, H, H0, intF, dF, W, and dW.
	[sigma,F,H,H0,intF,dF,W,dW] = trapF(options, bath);
	W(1)=1e100; %W(1) is the infinity, just make it huge, instead of the Inf, DJN 4/10/13
	W = W(:);

	%For no potential.
	%W=0*W;
	%dW=0*dW;


    dlambda= options.maxl/options.timesteps;
    dsigma = options.dsigma;
    alpha = options.bath.slope;
    g = options.g;

    n = length(sigma);

    dsigma2=dsigma*dsigma;   % Find dlambda^2 and dsigma^2
    dlambda2=dlambda*dlambda;

    A=sparse(n,n);

    A(1,1)=1;               % Define the matrix A, W and dW needed for our model
    for i=2:n-1
        A(i, i-1)=   -(    dlambda2/(dsigma2) - dlambda2/(2*dsigma)*W(i)                    );
        A(i, i)  = 1 -( -2*dlambda2/(dsigma2)                           + dlambda2*dW(i)    );
        A(i, i+1)=   -(    dlambda2/(dsigma2) + dlambda2/(2*dsigma)*W(i)                    );
    end

    A(n,n)=dsigma+dlambda;
    A(n,n-1)=-dlambda;

    %DJN
    %Define the initial profile and find Xmax.
    
    %Find the real Xmax.
    Max_H=interp1(sigma,H,options.maxsigma);
    xmax=Max_H/alpha;

    DJN_x=-[0:1:xmax];
    %DJN_eta=-0.0001/0.6065*exp(-2e-5*(1000+DJN_x).^2).*(1000+DJN_x); %alpha=0.01
    
    [DJN_eta DJN_flag]=eta_0(DJN_x);
    
    %DJN_eta=0.01*(1-tanh((1000+DJN_x)/200 ))/2
    
    %We need to convert (x, t, \eta, u) to (\sigma, \lambda, \phi, \psi)
    DJN_H=DJN_eta-DJN_x*alpha;
    DJN_Sigma=interp1(H, sigma, DJN_H);

    DJN_u=DJN_flag * DJN_eta.*sqrt(g./(-alpha*DJN_x));
    DJN_u(isnan(DJN_u))=0;







    DJN_Phi=2*g*DJN_eta;

    % Define the initial Phi (wave height)
    % Phi_nm1=-4*a*sigma.^(-1).*((sigma-s0)/p^2.*exp(-1*((sigma-s0)/p).^2)+(sigma+s0)/p^2.*exp(-((sigma+s0)/p).^2));
    % Phi_nm1(1)=0;

    % find when we need to switch from eta to linear theory boundry
    wasnon0=false;
    is0=false;
    time=0;
    counter=1;
    while(~(wasnon0&&is0))
        if(time>=40)
            if(eta_bound(Point)>0)
                counter=10^100;
                break
            else
            counter=-1;
            println('  - normal boundry')
            break
            end
        end
        Point=xmax+sqrt(abs(alpha*g*xmax))*time/(abs(alpha*g));
        if(~wasnon0)
            wasnon0=eta_bound(Point)>0;
        end
        is0=eta_bound(Point)==0;
        %if(wasnon0&&time>20)
         %   is0=eta_bound(Point)<10^-10;
        %end
            
        time=counter*dlambda;
        if ~(wasnon0&&is0)
            counter=counter+1;
        end
    end
    if(time<1.75)
        println('  - time to reach shore is too small. move away from shore')
        return
    end

    Phi_nm1=interp1(DJN_Sigma, DJN_Phi, sigma);
    Phi_nm1(isnan(Phi_nm1))=0;
    Phi_nm1(1)=0;

    Phi_nm1=Phi_nm1';   %Make it the column, DJN 4/10/13
    %Phi_sigma=Psi_lambda
    %Define the initial Psi and then the next time step (wave velocity)
    PSI_LAMBDA=zeros(n,1);                                                % Pre-allocate for speed

    PSI_LAMBDA(1)     =  (-Phi_nm1(3) + 4*Phi_nm1(2)-3*Phi_nm1(1))/(2*dsigma);     % Second order forwards difference
    PSI_LAMBDA(2:n-1) = (Phi_nm1(3:n) - Phi_nm1(1:n-2)) / (2*dsigma);           %phi(n+1)-Phi(n-1)
    PSI_LAMBDA(n)     =(-3*Phi_nm1(n) + 4*Phi_nm1(n-1)-Phi_nm1(n-2))/(2*dsigma); % Second order backwards difference



    u_sigma=interp1(DJN_Sigma, DJN_u, sigma);
    u_sigma(isnan(u_sigma))=0;
    
    Psi_nm1=(F.*u_sigma);
    Psi_nm1(isnan(Psi_nm1))=0;  % 1/0 is not good.
    Psi_nm1=Psi_nm1';   %Make it the column, DJN 4/10/13
    %zeros(n,1);                                          % psi=0, %Make it the column, DJN 4/10/13

    
    Psi_n=Psi_nm1+PSI_LAMBDA*dlambda;                                     % Compute psi at the second step
    %DJN 4/10/13 %Psi=Psi_n;                                                   % Define Psi as the nth step

    PHI_LAMBDA=zeros(n,1);
    % Phi_lambda=psi_sigma+W\psi
    %Find Phi at the next time step using Psi_n
    PHI_LAMBDA(1)=(-Psi_nm1(3)+4*Psi_nm1(2)-3*Psi_nm1(1))/(2*dsigma)+Psi_nm1(1)*W(1);     % Second order forwards difference
    PHI_LAMBDA(2:n-1) = (Psi_nm1(3:n)-Psi_nm1(1:n-2)) / (2*dsigma) + Psi_nm1(2:n-1) .* W(2:n-1); %psi(n+1)-psi(n+1)/(2dsigma)+psi(n)*W(n)
    PHI_LAMBDA(n)=(-3*Psi_nm1(n)+4*Psi_nm1(n-1)-Psi_nm1(n-2))/(2*dsigma)+Psi_nm1(n)*W(n); % Second order backwards difference
    
    Phi_n=Phi_nm1+dlambda*(PHI_LAMBDA);                                                   % Compute phi at the nth step
    %DJN 4/10/13 %Phi=Phi_n;                                                                   % Define Psi as the nth step



    %Pre-allocate for the speed %DJN 4/10/13
    n=length(sigma);
    dsigma = options.dsigma;

    Psiout=zeros(ceil(options.timesteps/options.keeprate), n);
    Phiout=zeros(ceil(options.timesteps/options.keeprate), n);
    lambda=zeros(ceil(options.timesteps/options.keeprate), 1);

    step=0;
    l=1;                            % Index to keep only parts of our informaion
    Phiout(l,:)=Phi_nm1;         % Keep the initial conditions
    Psiout(l,:)=Psi_nm1;
    lambda(l)=0;

    %DJN 4/10/13, we need to keep the second step too.
    step=1;
    if(mod(step,options.keeprate)==0) %Check if we need to keep it.
        l=2;                            % Index to keep only parts of our informaion
        Phiout(l,:)=Phi_n;           
        Psiout(l,:)=Psi_n;
        lambda(l)=step*dlambda;
    end

    figure(1); clf;
    for step=2:options.timesteps    %we start from the third step, since the first two are already computed, DJN 4/10/13
        %DJN  b(1)=0;                     % Define b as the right side of our system
        %     for i=2:n-1
        %         b(i)=2*Psi_n(i)-Psi_nm1(i);
        %     end
        %     b(n)=0;
        
        b=2*Psi_n-Psi_nm1;          %Convert into the vector operation, DJN 4/10/13
        b(1)=0;
        
        if(counter>step&&counter~=-1)% if we have a moving boundry
            %set phi equal to F U where u=\sqrt(g/h)eta and h=alpha*xmax
            A(n,n)=1;
            A(n,n-1)=0;
            %eta at the x-t point(assuming linear velosity
            b(n)=interp1(sigma,F,maxsigma) * sqrt(g/(alpha*xmax)) * eta_bound(xmax+sqrt(abs(alpha*g*xmax))*step*dlambda/(abs(alpha*g)));
        else
            %implicit method
            % we have psi_lambda=-psi_sigma
            A(n,n)=dsigma+dlambda;
            A(n,n-1)=-dlambda;
            b(n)=dsigma*Psi_n(n);
        end

        %DJN  Psi=A\b;                    %solve for Psi and set the timesteps up one
        %     Psi_nm1=Psi_n;
        %     Psi_n=Psi;
        Psi_nm1=Psi_n;              %We don't really need Psi vector, it just got eliminated to save time, DJN 4/10/13
        Psi_n=A\b;
        
        if(counter>step&&counter~=-1)% if we have a moving boundry
            Psi_n(end)=2*g*eta_bound(xmax+sqrt(abs(alpha*g*xmax))*step*dlambda/(abs(alpha*g)));
        end
        
        
        
        PHI_LAMBDA(1)     = (  -Psi_n(3) + 4*Psi_n(2)-3*Psi_n(1))/(2*dsigma)+W(1)*Psi_n(1);     % Second order forwards difference
        PHI_LAMBDA(2:n-1) = ( Psi_n(3:n) - Psi_n(1:n-2)) ./ (2*dsigma) + W(2:n-1).*Psi_n(2:n-1);  %psi(n+1)-psi(n+1)/(2dsigma)+psi(n)*W(n)
        PHI_LAMBDA(n)     = (-3*Psi_n(n) + 4*Psi_n(n-1)-Psi_n(n-2))/(2*dsigma)+W(n)*Psi_n(n); % Second order backwards difference

        Phi=4/3*Phi_n-1/3*Phi_nm1+2/3*PHI_LAMBDA*dlambda;                             % Define the next Phi
        Phi_nm1=Phi_n;
        Phi_n=Phi;
        
        if(mod(step,options.keeprate)==0)              % Keep information at some points
            Psiout(l,:)=Psi_n;              % save the values at the current time step (written into the *_n arrays)
            Phiout(l,:)=Phi_n;
            lambda(l)=step*dlambda;
            
            plot(sigma(1,1:n-2),Phiout(l,1:n-2),'.b')
            set(gca,'xdir','reverse')
            title(['Step ' num2str(step) ' (' num2str(100 *step /options.timesteps) '%)'])
            drawnow();
            l=l+1;
        end
    end


	results = struct('phi',Phiout, 'psi',Psiout, 'lambda',lambda ,'x0',DJN_x ,'eta0',DJN_eta, 'F',F, 'intF',intF);
end %function
