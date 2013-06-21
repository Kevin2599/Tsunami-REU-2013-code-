function y=eta_bound(x)
	a=.5;
	x_l=7000;
	l=5;
	y=a*exp(-((x-x_l)/l).^2);% gause wave
end