%[sigma,F,H,H0,intF,dF,W,dW] = trapF(beta,y0,dsigma,sigmamax,ymax,g=1);
%
%Given a symmetric trapezoidal cross-section defined by beta and y0, outputs the
%variety of bayometrically determined variables required by our model.
%
%Inputs:
%-beta is the slope of the trapezoid walls.
%-y0 is half the length of the trapezoid's bottom.
%-dsigma is the desired interval of the output sigma vector. Should not be
%smaller than 0.001.
%-sigmamax is the desired maximum of the output sigma vector.
%-ymax is a parameter that must be sufficiently large to allow for the
%desired sigmamax.
%-g is the value assigned to gravity. The default is 1.
%
%Outputs:
%-sigma is the independent output variable, with step size dsigma and
%maximum sigmamax.
%-F(sigma) is the smoothed desired function. It is calculated internally
%with a sigma step size of 0.001, smoothed, then, along with the remaining
%outputs, is cut down to use the sigma step size specified by dsigma. F is
%the stitched composite of two functions: close to the shore, it uses an
%analytic asymptotic expression for F, while farther out to sea it uses
%numeric methods. The stitching point is dynamically determined to be the
%point at which those two functions intersect.
%-H(sigma) is the total depth of water in the center of the channel, as a
%function of sigma.
%-H0(sigma) is the effective depth of water in the channel. It is
%calculated by dividing the area of the water in a cross-section by the
%width of the channel at the surface.
%-intF(sigma) is the integration of F(sigma). This is calculated using the
%identity intF = 2*g*H(sigma).
%-dF(sigma) is the smoothed first derivative of F(sigma). It is calculated using a
%second order central difference scheme, then smoothed using the matlab
%smooth() function. There is a jump at the stitching point between the
%asymptotic and numeric functions, which is smoothed using the fixit()
%function, which is a series of repeated linear smoothings. Thus, running
%trapF requires that the fixit.m function is accessible.
%-W(sigma) is a function of F and dF.
%-dW(sigma) is the derivative of W with respect to sigma.
%
%If you have questions, you can contact me (Lander Ver Hoef) at
%lverhoef@alaska.edu.

function [sigma,F,H,H0,intF,dF,W,dW] = trapF(options,ymax)

	beta     = options.bath.trap_slope;
	y0       = options.bath.trap_width/2;
	dsigma   = options.dsigma;
	sigmamax = options.maxsigma;
	g        = options.g;


	%We define our internal intervals.
	yinterval=.001;
	sigmainterval=0.001;

	%We do only half of y since it's symmetric.
	y = 0:yinterval:ymax; 
	%Next let's find the index of y0.
	y0index = find(y==y0,1);
	%We pre-allocate f, then use a loop to define the piecewise trapezoidal
	%cross-section.
	f=zeros(1,length(y));
	f(y0index:end) = beta*( y(y0index:end) - y0 );


	%We need to find h0temp.  This is defined as:
	%S/(dS/dH) = S/(ap - an).  Thus, for our existing vectors y and f (where f
	%acts like H) we need to find S and ap and an.  Since we are symmetric,
	%an=-ap.
	%For a given value of H (or rather, f) the value for ap is y.

	%To find S, we need to integrate.
	%S = H*(ap-an)- int_an^ap (f) \dy
	%That integral is just 2*int_0^ap f dy

	%We numerically perform this integral with Simpson's rule.
	%{
	for i = 1:length(f);
	    S(i)=f(i)*2.*y(i)-2*quad(f,0,y(i));
	end
	%}
	S=f*2.*y-2*cumsimps(y,f);

	%We introduce a loop that removes values from horizontal planes in f(y) so
	%that it is invertible.

	%%DJN
	% % i=2;
	% % while i<=length(y)
	% %     if f(i)==f(i-1)
	% %         f(i-1)=[];
	% %         y(i-1)=[];
	% %         S(i-1)=[];
	% %     else
	% %         i=i+1;
	% %     end
	% % end

	[DJN_f,DJN_m,DJN_n]=unique(f);
	f=f(DJN_m);
	y=y(DJN_m);
	S=S(DJN_m);


	%Now, we find h0temp using the definition h0 = S / (ap - an)
	h0temp = S ./ (2*y);

	%Next: find our temp sigma as a function of H and h0, here represented by f
	%and h0temp.

	integrand = sqrt(g./h0temp);
	integrand(1) = 0;
	sigmatemp = cumsimps(f, integrand);
	%{
	sigmatemp=zeros(1,length(f));
	for i = 1:length(f);
	    sigmatemp(i) = quadv(integrand,0,f(i));
	end
	%}

	%Now, we define our longer lasting variables! We use a loop to check if the
	%sigmamax that we want is within the acceptable interval.
	if sigmamax<max(sigmatemp)
	    preciseSigma = 0:sigmainterval:sigmamax;
	else
	    preciseSigma = 0:sigmainterval:max(sigmatemp);
	    disp('Warning: The maximum of sigmatemp is smaller than the specified sigmamax. Increase ymax.')
	end

	%We interpolate values for H as a function of preciseSigma based on our existing
	%data from sigmatemp and f.
	preciseH = interp1(sigmatemp, f, preciseSigma);

	%Now, we find H0 based on our existing data from h0temp and f.
	preciseH0 = interp1(f, h0temp, preciseH);

	%Now that we have H0(preciseSigma), we can easily calculate the rough F!
	roughF = 2*sqrt(g*preciseH0);

	%Now let's introduce the asymptotics.
	%For the trapezoidal case, we know that F ~ sigma -
	%(1/(12*g*beta*y0))*sigma.^3 + (1/(48*g*beta^2*y0^2))*sigma.^5 as sigma
	%goes to 0.

	asympF = preciseSigma -(1/(12*g*beta*y0))*preciseSigma.^3 + (1/(48*g^2*beta^2*y0^2))*preciseSigma.^5;

	%Now we need to stitch these together.

	%Let's smooth the numeric F
	numericF = smooth(roughF,15)';

	%Because we have smoothed F, we can find the point at which the asymptotic F equals the numeric F
	stitchPoint = interp1(numericF(20:end) - asympF(20:end),preciseSigma(20:end),0);
	%Find stitching point in terms of the indices
	for i=2:length(preciseSigma)
	    if preciseSigma(i) > stitchPoint
	        stitchIndex = i;
	        break;
	    end
	end

	%We do a piecewise stitching to get an unsmoothed F.
	preciseF = zeros(1,length(numericF));
	preciseF(1:(stitchIndex-1)) = asympF(1:(stitchIndex-1));
	preciseF(stitchIndex:length(numericF)) = numericF(stitchIndex:length(numericF));

	%Now to find the integral of F using the relation int F = 2gH.
	tempIntF=2*g*preciseH;

	%We change preciseSigma into sigma with step size dsigma for the larger
	%modeling program.
	sigma=0:dsigma:max(preciseSigma);
	H = zeros(1,length(sigma));
	H0 = zeros(1,length(sigma));
	intF = zeros(1,length(sigma));
	j=1;

	%Now we go from our precise variables with a step size of 0.001 in sigma to
	%the output variables with step size dsigma.
	for i = 1:length(preciseSigma)
	    if ~mod(i-1,dsigma/sigmainterval)
	       sigma(j) = preciseSigma(i);
	       F(j) = preciseF(i);
	       H(j) = preciseH(i);
	       H0(j) = preciseH0(i);
	       intF(j) = tempIntF(i);
	       j = j+1;
	    end
	end

	%And now we find the derivative of F with respect to sigma
	%First, we pre-allocate dF and define dF(1) using a forwards difference
	%scheme
	n=length(sigma);
	dF=zeros(1,n);
	dF(1)=(-F(3)+4*F(2)-3*F(1))/(2*dsigma);

	%Use a second order central difference scheme 
	dF(2:n-1) = (F(3:n)-F(1:n-2)) / (2*dsigma);

	%Find the endpoint using a backwards difference scheme
	dF(n)=(-3*F(n)+4*F(n-1)-F(n-2))/(2*dsigma);

	%Now we find the index that corresponds to the stitching point in the new
	%vectors.
	sigStitchIndex = find(sigma==preciseSigma(stitchIndex),1);
	%Now let's smooth our dF
	smoothdF = (fixit(smooth(dF,15),sigStitchIndex-40,sigStitchIndex+40,30))';

	%corect last few points
	for temp=2:length(smoothdF)
	   if abs(smoothdF(temp-1)-smoothdF(temp))>.03
	       smoothdF(temp)=smoothdF(temp-1);
	   end
	end

	%Now find W using F and dF.
	W=(2-smoothdF)./F;



	%plot(smoothdF(1:length(smoothdF)-1)-smoothdF(2:length(smoothdF)),'.k')





	%Now let's define our asmyptotic W.
	asympW = 1./sigma + sigma./(3*g*beta*y0);

	%Find dW using 
	dW=-(1./F).^2.*(2*ones(1,length(dF))-smoothdF);
	%Let's also do the asmyptotics of dW.
	asympdW = -(1./(sigma.^2)) + 1/(3*g*beta*y0);
	dF = smoothdF;
end
