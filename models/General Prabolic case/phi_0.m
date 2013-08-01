function [phi pot,A, p, sigma_0] =phi_0(sigma)
% 	A=.1;
%     p=1;
%     sigma_0=2.9;
    A=.5;
    p=1.5;
    sigma_0=15;
	phi=-4*A*sigma.^(-1).*((sigma-sigma_0)/p^2.*exp(-1*((sigma-sigma_0)/p).^2)+(sigma+sigma_0)/p^2.*exp(-((sigma+sigma_0)/p).^2));
    phi(abs(phi)<1e-5)=0;
    phi(1)=0;
	pot=0;
    plot(sigma,phi,'.b')
end