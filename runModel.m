function [Phiout Psiout lambda] = runModel(timesteps,keeprate,sigma,Phi_nm1,Phi_n,Psi_nm1,Psi_n,counter,A,dlambda,dsigma,W,PHI_LAMBDA,F)
    %Pre-allocate for the speed %DJN 4/10/13
    n=length(sigma);

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

    figure(1); clf;
    for step=2:timesteps    %we start from the third step, since the first two are already computed, DJN 4/10/13
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
        
        if(mod(step,keeprate)==0)              % Keep information at some points
            Psiout(l,:)=Psi_n;              % save the values at the current time step (written into the *_n arrays)
            Phiout(l,:)=Phi_n;
            lambda(l)=step*dlambda;
            
            plot(sigma(1,1:n-2),Phiout(l,1:n-2),'.b')
            title(['Step ' num2str(step) ' (' num2str(100 *step /timesteps) '%)'])
            drawnow();
            l=l+1;
        end
    end
end