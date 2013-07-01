function y=eta_bound(x)
	a=.0005;
	x_l=8000;
	l=2000;
	y=a*exp(-((x-x_l)/l).^2);% gause wave
    y(y<10^-5)=0;
	y(y~=0)=0;
   % y=.2/14.11;
end