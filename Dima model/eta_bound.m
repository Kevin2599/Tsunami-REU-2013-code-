function [eta_boundSigma on]=eta_bound(sigma)
a=.0005;
sigma_l=9000;
l=2000;
eta_boundSigma=a*exp(-((sigma+sigma_l)/l).^2);% gause wave
eta_boundSigma(eta_boundSigma<10^-5)=0;
eta_boundSigma(eta_boundSigma~=0)=0;
% eta_boundSigma=.2/14.11;
on=1;
if sum(abs(eta_boundSigma))==0
    on=0;
end
end