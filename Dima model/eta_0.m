function [y U_on] =eta_0(x)


y=-9.0315e-4*exp(-1.5e-5*(1000+x).^2).*(1000+x);% N-wave \alpha =0.05
 
    %  y=0.25 + 0.25*tanh((-x - 1000)/15);
    %y=.5./(1+exp(x+1000));
   % y=.5*normcdf(-x,1000,250);
	y(abs(y)<1e-5)=0;
	% y=y*0;
	U_on = 0;
end