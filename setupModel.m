function [Phi_nm1,Phi_n, Psi_nm1,Psi_n,counter,A,dlambda,W,PHI_LAMBDA, DJN_x,DJN_eta] = setupModel(sigma,W,dW,H,F,options)
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

end