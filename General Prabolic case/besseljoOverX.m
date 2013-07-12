function y=besseljoOverX(a,lambda,x)
% k=0:50;
%  (-1).^k(end)./(gamma(k(end)+a+1).*factorial(k(end))).*(1/2).^(2*k(end)+a).*(sqrt(-lambda)).^(2*k(end)+a).*x.^(2*k(end));
% y=sum((-1).^k./(gamma(k+a+1).*factorial(k)).*(1/2).^(2*k+a).*x.^(2*k).*(sqrt(-lambda)).^(2*k+a));
if x==0
    y=(2^(-1-a)*(sqrt(-lambda))^(a+1))/(gamma(2+a));

else
    y=NaN;
end
end