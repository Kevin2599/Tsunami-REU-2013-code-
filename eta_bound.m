function y=eta_bound(x)
	a=.05;
	x_l=8000;
	l=100;
	y=a*exp(-((x-x_l)/l).^2);% gause wave
	%y(:,:)=0;
end