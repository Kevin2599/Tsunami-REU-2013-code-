% This program will find the solution to the tsunami runup problem on a
% parabolic beach (cross sections have form |y|^m) with a constant slope in
% the x direction:
% The following varables are used in this program:
% W        - Vector that is used to find A.
% A        - n by n matrix that is used to solve the wave equation.
% I        - Matrix to look for breaks in time.
% a        - The amplitud of our gauss pulse.
% alpha    - The slope of th beach.
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
% m        - |y|^m.
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
% Psi_n    - Length n vector that holds the curent time step for our solution
%            to the wave equation. Needed to shuffle data.
% Psi_nm1  - Length n vector that holds the curent time step for our solution
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
% timestpes- Sets the change in lambda.
% u1       - n by length(lambda) matrix for our velocity output fot the
%            exact solution.
% u2       - n by length(lambda) matrix for our velocity output fot the
%            aprox solution.
% x1       - Our distance output for the exact solution. NOTE MATRIX
% x2       - Our distance output for the aprox solution. NOTE MATRIX

function results = evalUnboundedParabolic(options)

    m  = 2;                     % m difines the bay shape |y|^m
    a  = .5;                    % a is the amplutude of our pulse
    s0 = 15;                   % so is the mean of out pulse
    p  = 1.5;                   % p is the  varence in pulse

    maxl      = options.maxl
    timesteps = options.timesteps;            % number of time steps between \lambda=0, and \lambda=maxl, %DJN 4/10/13
    keeprate  = options.keeprate;     % keep every \it{keeprate}-th step.
    plotb     = 1;                 % Bool to plot
    dlambda   = options.maxl/options.timesteps;
    dsigma    = options.dsigma;
    alpha     = options.bath.slope;
    g         = options.g;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build starting need information from user inputs and build the matrix A
    % that will be used to solve our system

    disp('Building model...')



    dsigma2=dsigma*dsigma;   % Find dlambda^2 and dsigma^2
    dlambda2=dlambda*dlambda;
    sigma=0:dsigma:options.maxsigma;       % Define sigma
    n=length(sigma);

    A=sparse(n,n);
    b=zeros(n,1);


    A(1,1)=1;               % Define the matrix A, W and dW needed for our model
    for i=2:n-1
        W(i)  = (m+2)/(m*sigma(i));
        dW(i) = -(m+2)/(m*sigma(i)^2);
        
        A(i, i-1)=   -(    dlambda2/(dsigma2) - dlambda2/(2*dsigma)*W(i)                    );
        A(i, i)  = 1 -( -2*dlambda2/(dsigma2)                           + dlambda2*dW(i)    );
        A(i, i+1)=   -(    dlambda2/(dsigma2) + dlambda2/(2*dsigma)*W(i)                    );
    end
    W(n)  = (m+2)/(m*sigma(n));
    dW(n) = -(m+2)/(m*sigma(n)^2);
    A(n,n)=1;



    % define the initial Phi (wave height)
    Phi_nm1=-4*a*sigma.^(-1).*((sigma-s0)/p^2.*exp(-1*((sigma-s0)/p).^2)+(sigma+s0)/p^2.*exp(-((sigma+s0)/p).^2)); % at \lambda=0
    %Phi_nm1=a*(exp(-1*((sigma-s0)/p).^2));
    %Phi_nm1=2*g*-0.0001/0.6065*exp(-2e-5*(1000+x).^2).*(1000+x);
    Phi_nm1(1)=0;



    %Define the initial Psi and then the next time step (wave velosity)
    G=zeros(1,n);                                                % Prealacate for speed
    G(1)=(-Phi_nm1(3)+4*Phi_nm1(2)-3*Phi_nm1(1))/(2*dsigma);     % Second order forwards differene
    for i=2:n-1
        G(i)=(Phi_nm1(i+1)-Phi_nm1(i-1))/(2*dsigma);             % Second order centeral difference
    end
    G(n)=(-3*Phi_nm1(n)+4*Phi_nm1(n-1)-Phi_nm1(n-2))/(2*dsigma); % Second order backwards difference

    Psi_nm1=zeros(1,n);                                          % psi=0  at \lambda=0


    Psi_n=Psi_nm1+G*dlambda;                                     % Compute psi at the second step, \lambda=d\lambda, It is the second order accurate since \psi(0)=0
    Psi=Psi_n;                                                   % Define Psi as the nth step


    %find Phi at the next time step using Psi_n
    G(1)=(-Psi_nm1(3)+4*Psi_nm1(2)-3*Psi_nm1(1))/(2*dsigma)+Psi_nm1(1)*W(1);     % Second order forwards difference
    for i=2:n-1
        G(i)=(Psi_nm1(i+1)-Psi_nm1(i-1))/(2*dsigma)+Psi_nm1(i)*W(i);             % Second order centeral difference
    end
    G(n)=(-3*Psi_nm1(n)+4*Psi_nm1(n-1)-Psi_nm1(n-2))/(2*dsigma)+Psi_nm1(n)*W(n); % Second order backwards difference
    Phi_n=Phi_nm1+dlambda*(G);                                                   % Compute phi at the nth step
    Phi=Phi_n;                                                                   % Define Psi as the nth step

    %\psi_\lambda=\phi_\sigma
    %\phi_\lambda=\psi_\sigma+W*\psi
    %\phi_\lambda_\lambda=\psi_\lambda_\sigma+W*\psi_\lambda=   \phi_\sigma_\sigma+W*\phi_\sigma
    %\phi(dL)=\phi(0)+dL*\frac{\phi}{\lambda}(0)+
    %0.5*dL^2*\frac{\phi^2}{\lambda^2}(0)    %the last term is missing, but it


















    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the model for Psi and Phi

    disp('Running model...')
    l=1;                            % Index to keep only parts of our informaion
    Psiout(:,l)=Psi_nm1(:);         % Keep initial conditions
    Phiout(:,l)=Phi_nm1(:);
    lambda(l)=0;
    l=l+1;
    for step=1:maxl*timesteps
        
        b(1)=0;                     % Define b as the right side of our system
        for i=2:n-1
            b(i)=2*Psi_n(i)-Psi_nm1(i);
        end
        b(n)=0;
        
        Psi=A\b;                    %solve for Psi and set the timesteps up one
        Psi_nm1=Psi_n;
        Psi_n=Psi;
        
        G(1)=(-Psi_n(3)+4*Psi_n(2)-3*Psi_n(1))/(2*dsigma)+W(1)*Psi_n(1);     % Second order forwards differene
        for i=2:n-1
            G(i)=(Psi_n(i+1)-Psi_n(i-1))/(2*dsigma)+W(i)*Psi_n(i);           % Second order centeral differene
        end
        G(n)=(-3*Psi_n(n)+4*Psi_n(n-1)-Psi_n(n-2))/(2*dsigma)+W(n)*Psi_n(n); % Second order backwards differene
        Phi=4/3*Phi_n-1/3*Phi_nm1+2/3*G*dlambda;                             % Define the next Phi
        Phi_nm1=Phi_n;
        Phi_n=Phi;
        
        
        if(step*dlambda*keeprate==floor(step*dlambda*keeprate))                          % Keep information at some points
            Psiout(:,l)=Psi(:);
            Phiout(:,l)=Phi(:);
            lambda(l)=step*dlambda;
            l=l+1;
        end
    end
























    clearvars -except 'Phiout' 'Psiout' 'a' 'p' 's0' 'sigma' 'lambda' 'm' 'g' 'alpha' 'plotb'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the exact solutions
    %WARNING: If initial Psi is non zero then the exact solution is not our
    %solution. Also if m~=2 then the exact solution is not our solution this
    %program will automatica
    Exact=0;
    if ((m==2)&&(~max(abs(Psiout(:,1)))))
        Exact=1;
    end

    if Exact
        disp('Finding Explicit Data...')
        % Define a mesh with the corect lambda for exact solution
        lambda=-lambda;
        [LAM, SIG] = meshgrid(lambda, sigma);
        
        % Find exact phi and psi.
        phi = a./ SIG .* ( -2.*(SIG+LAM-s0)/(p^2).* exp(-(((SIG+LAM-s0)/p).^2)) -2.*(SIG-LAM-s0)/(p^2).* exp(-(((SIG-LAM-s0)/p).^2)) -2.*(SIG+LAM+s0)/(p^2).* exp(-(((SIG+LAM+s0)/p).^2))-2.*(SIG-LAM+s0)/(p^2).* exp(-(((SIG-LAM+s0)/p).^2)));
        psi = -a./ SIG .* ( -(1./SIG + 2*(SIG + LAM - s0)/(p^2)).* exp(-(((SIG + LAM - s0)/p).^2)) + (1./SIG + 2*(SIG - LAM - s0)/(p^2)).* exp(-(((SIG - LAM - s0)/p).^2)) - (1./SIG + 2*(SIG + LAM + s0)/(p^2)).* exp(-(((SIG + LAM + s0)/p).^2)) + (1./SIG + 2*(SIG - LAM + s0)/(p^2)).* exp(-(((SIG - LAM + s0)/p).^2)));
      


    end
    if ~Exact
        disp('No Explicit Data...')
    end



















    clear LAM SIG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert back to physical varables

    if Exact
        disp('Converting Exact and Approx data...')
        %rescaile phi and psi for both cases to convert
        phi=m/(m+1)*phi;
        psi=m/(m+1)*psi;
        Phiout=m/(m+1)*Phiout;
        Psiout=m/(m+1)*Psiout;
        
        
        % Define the bayomatry F and intF
        F = (m/(m+1))*sigma;
        intF = m/2*(1/(m+1)*(sigma.^2));
        
        % Data Needed to convert both exact and aprox data
        [LAM, Fgrid] = meshgrid(-lambda, F);
        [dummy, intgrid] = meshgrid(-lambda, intF);
        clear dummy
        
        % Convert Exact
        u1 = psi./Fgrid;
        eta1=(phi-u1.^(2))/(2*g);
        t1=((LAM-u1)/(alpha.*g));
        x1 = (phi-u1.^(2)-intgrid)/(2*alpha*g);
        t1=abs(t1);
        
        % Convert Aprox.
        u2 = Psiout./Fgrid;
        eta2=(Phiout-u2.^(2))/(2*g);
        t2=((LAM-u2)/(alpha.*g));
        x2 = (Phiout-u2.^(2)-intgrid)/(2*alpha*g);
        t2=abs(t2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clear all data that is nolonger needed
        clearvars -except  'u1' 'u2' 'eta1' 'eta2' 'x1' 'x2' 't1' 't2' 'm' 'g' 'alpha' 'lambda' 'Exact' 'plotb'
    end


    if ~Exact
        disp('Converting Approx data...')
        %rescaile phi and psi for both cases to convert
        Phiout=m/(m+1)*Phiout;
        Psiout=m/(m+1)*Psiout;
        
        % Define the bayomatry F and intF
        F = (m/(m+1))*sigma;
        intF = m/2*(1/(m+1)*(sigma.^2));
        
        % Data Needed to convert both exact and aprox data
        [LAM, Fgrid] = meshgrid(-lambda, F);
        [dummy, intgrid] = meshgrid(-lambda, intF);
        clear dummy
        
        % Convert Aprox.
        u2 = Psiout./Fgrid;
        eta2=(Phiout-u2.^(2))/(2*g);
        t2=((LAM-u2)/(alpha.*g));
        x2 = (Phiout-u2.^(2)-intgrid)/(2*alpha*g);
        t2=abs(t2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clear all data that is nolonger needed
        clearvars -except  'u2'  'eta2' 'x2' 't2' 'm' 'g' 'alpha' 'lambda' 'Exact' 'plotb'
    end









    %{
    format long
     [maxi,W1]=max(eta1(2,:))
      [mini,W2]=min(eta1(2,:))
      [maxi,W11]=max(eta2(2,:));
      [mini,W21]=min(eta2(2,:));
    x=[-100,100];

    plot(x1(:,W1),eta1(:,W1), 'r')
    hold on
    plot(x1(:,W11),eta2(:,W11),'c')
    plot(x,alpha*x)
    plot(0,0,'^b')




    plot(x1(:,W2),eta1(:,W2), 'r')
    plot(x1(:,W21),eta2(:,W21),'c')
    plot(x,alpha*x)
    plot(0,0,'^b')
    axis((1+floor(max(max(abs(eta2)))))*[-10*max(max(x2)) 1.5*max(max(x2)) alpha*-10*max(max(x2)) 4*alpha*max(max(x2))])
    leg=legend('Exact solution', 'modeled solution','Bay','Initial water level','Exact solution', 'modeled solution','Location','Best');
    xlabel('x')
    return;
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the data

    % Look for break in time.
    disp('Plotting...')
    toc
    [dummy,I]=sort(t2*alpha,2);
    found=0;
    brokeat=length(t2(1,:))+1;
    for j=1:length(t2(1,:))
        if I(2,j)~=j
            found=1;
            brokeat=j;
            break
        end
    end


    if plotb
        if Exact
            % Plot to look for global error and information
            slope=zeros(1,length(lambda));
            breakc=slope;
            for j=1:length(lambda)
                for i=1:length(eta1(:,1))-1
                    slope(i)=(eta1(i+1,j)-eta1(i,j))/(x1(i+1,j)-x1(i,j));
                end
                breakc(j)=max(slope(:));
            end
            for i=1:length(lambda)
                if ((breakc(i)>=1/2*alpha)||(i==brokeat))
                    disp('BROKE...')
                    if found
                        disp('Numerical')
                    end
                    break
                end
                plot(eta1(:,i), 'r')
                hold on
                plot(abs(eta1(:,i)-eta2(:,i)),'b')
                hold off
                axis([0 300 min(min(eta1)) max(max(eta1))])
                leg=legend('Exact solution', 'Error');
                %set(leg,'Location','Best')
                xlabel('Sigma')
                title(num2str(t1(2,i)))
                pause(0.01)
            end
            
            % Plot at the shore
            x=[-100,100];
            for i=1:5:length(lambda)
                if ((breakc(i)>=1/2*alpha)||(i==brokeat))
                    disp('BROKE...')
                    if found
                        disp('Numerical')
                    end
                    break
                end
                plot(x1(:,i),eta1(:,i), 'r')
                hold on
                plot(x1(:,i),abs(eta1(:,i)-eta2(:,i)),'c')
                plot(x,alpha*x)
                plot(0,0,'^b')
                hold off
                axis((1+floor(max(max(abs(eta2)))))*[-10*max(max(x2)) 1.5*max(max(x2)) alpha*-10*max(max(x2)) 4*alpha*max(max(x2))])
                leg=legend('Exact', 'Error','Location','Best');
                %set(leg)
                xlabel('x')
                title(num2str(t1(2,i)))
                pause(0.01)
            end
        end
        
        
        if ~Exact
            % Plot to look for global error and information
            slope=zeros(1,length(lambda));
            breakc=slope;
            for j=1:length(lambda)
                for i=1:length(eta2(:,1))-1
                    slope(i)=(eta2(i+1,j)-eta2(i,j))/(x2(i+1,j)-x2(i,j));
                end
                breakc(j)=max(slope(:));
            end
            for i=1:length(lambda)
                if ((breakc(i)>=1/2*alpha)||(i==brokeat))
                    disp('BROKE...')
                    if found
                        disp('Numerical')
                    end
                    break
                end
                plot(eta2(:,i), 'r')
                axis([0 300 min(min(eta2)) max(max(eta2))])
                leg=legend('Aprox solution');
                % set(leg,'Location','Best');
                xlabel('Sigma')
                title(num2str(t2(2,i)))
                pause(0.01)
                
            end
            
            % Plot at the shore
            x=[-100,100];
            for i=1:length(lambda)
                if ((breakc(i)>=1/2*alpha)||(i==brokeat))
                    disp('BROKE...')
                    if found
                        disp('Numerical')
                    end
                    break
                end
                plot(x2(:,i),eta2(:,i), 'r')
                hold on
                plot(x,alpha*x)
                plot(0,0,'^b')
                hold off
                axis((1+floor(max(max(abs(eta2)))))*[-10*max(max(x2)) max(max(x2)) alpha*-10*max(max(x2)) alpha*max(max(x2))])
                leg=legend('Aprox Solution');
                %set(leg,'Location','Best')
                xlabel('x')
                title(num2str(t2(2,i)))
                pause(0.01)
            end
        end
    end