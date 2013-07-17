function [y oneorzeroU oneorzeroEta] =eta_0(x)
	%y=-9.0315e-4*exp(-1.5e-5*(1000+x).^2).*(1000+x);% gause wave
	%y(abs(y)<1e-3)=0;
	%y(y~=0)=0;


	y=-0.0001/0.6065*exp(-2e-5*(1000+x).^2).*(1000+x); 
	%y=.1*exp(-((x+1000)/50).^2);
    
	y(abs(y)<1e-5)=0;
	% y=y*0;
	oneorzeroU = 0;
    oneorzeroEta=1;
end