function [Phi,x_lin t_lin eta_lin u_lin]=parasolution()
Num_Roots=100;
m=2;

lambda_S=0;
lambda_E=50;
%plotcoef=0;


alpha=.05;
x=-(0.000001:.05:400);
g=9.81;
sigma=sqrt(-x*alpha)*2*g*(m+1)/m;
Sigma_L=max(sigma);

lambda=(lambda_S:lambda_E/200:lambda_E);

disp('bessel')
Zeros_of_bessel=besselzero((1/m)+1,Num_Roots-1,1);
Lambda_n=-([0 Zeros_of_bessel']./Sigma_L).^2;% build the eigenvalues from zeros of bessel
if(m==2)
    plot(tan(Sigma_L*sqrt(-Lambda_n))-Sigma_L*sqrt(-Lambda_n))%compares eigen values to expected ones for m=2, shows deviation from 0 for each eigen value.
    pause(1)
end

%find constants
disp('consts')
C_n=zeros(Num_Roots,1);
D_n=zeros(Num_Roots,1);
for n=1:Num_Roots%need to make besselj better so that x can be 0
    
    besselpart=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
    
    plot(besselpart)
    pause(.1)
    
    C_n(n)=trapz(sigma,(sqrt(-1*Lambda_n(n))*2*g*eta_0(-sigma).*besselpart));
    D_n(n)=0.*trapz(sigma,(besselpart.*cumsimps(sigma,m/(m+1).*sigma.*eta_0(-sigma).*sqrt(g./(alpha*-x)))));
    %D_n(n)=0;
    
end

%if plotcoef
%    plot(C_n)
%    hold on
%    plot(D_n,'r')
%    return;
%end

NEED_ROOTS=max([abs(max(C_n(floor(Num_Roots*.9):end))),abs(min(C_n(floor(Num_Roots*.9):end))),abs(max(D_n(floor(Num_Roots*.9):end))),abs(min(D_n(floor(Num_Roots*.9):end)))]);%measure of how far from 0 last 10% of constants are

if NEED_ROOTS>10
    disp('Need More Roots')
    plot(C_n)
    hold on
    plot(D_n,'r')
    hold off
    
    while NEED_ROOTS>10 %adds more roots so that eigen expansion is more accurate if need be; make more restrictive for larger values of xmax.
        Num_Roots=Num_Roots+1;
        temp=besselzero((1/m)+1,Num_Roots-1,1);
        Lambda_n(Num_Roots)=-1*(temp(end)./Sigma_L).^2;%add the new lambda
        n=Num_Roots;
        besselpart=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
        C_n(n)=trapz(sigma,(sqrt(-1*Lambda_n(n))*2*g*eta_0(-sigma).*besselpart));
        D_n(n)=0.*trapz(sigma,(besselpart.*cumsimps(sigma,m/(m+1).*sigma.*eta_0(-sigma).*sqrt(g./(alpha*-x)))));
        %D_n(n)=0;
        NEED_ROOTS=max(max(C_n(floor(Num_Roots*.9):end))-min(C_n(floor(Num_Roots*.9):end)),max(D_n(floor(Num_Roots*.9):end))-min(D_n(floor(Num_Roots*.9):end)));
        plot(C_n)%visual aid to see last portion of constants approaching 0 as new roots are added
        hold on
        plot(D_n,'r')
        hold off
        pause(.0000000000000000001)
    end
    
    if(m==2)%compares eigen values to expected ones for m=2, shows deviation from 0 for each eigen value
        plot(tan(Sigma_L*sqrt(-Lambda_n))-Sigma_L*sqrt(-Lambda_n))
        pause(1)
    end
    
end



disp('put together')
Phi= zeros(Num_Roots,length(sigma),length(lambda));
for j=1:length(lambda)
    besselpart=besselj(1/m,sqrt(-1*Lambda_n(1))*sigma)./sigma.^(1/m);
    Phi(1,:,j)=besselpart.*(C_n(1).*sin(sqrt(-1*Lambda_n(1))*lambda(j))+D_n(1).*cos(sqrt(-1*Lambda_n(1))*lambda(j)));%both D_n had sin instead of cos YIKES, look through paper work (finally)
end
for n=2:Num_Roots
    besselpart=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
    for j=1:length(lambda)
        Phi(n,:,j)=Phi(n-1,:,j)+besselpart.*(C_n(n).*sin(sqrt(-1*Lambda_n(n))*lambda(j))+D_n(n).*cos(sqrt(-1*Lambda_n(n))*lambda(j)));
    end
    disp(Num_Roots-n)
end



disp('finding PHI&Psi')

p2(:,:)=Phi(Num_Roots,:,:);
phi=zeros(length(sigma),length(lambda)-1);
psi=zeros(length(sigma)-1,length(lambda));

phi(:,1)=(p2(:,2)-p2(:,1))/(lambda(2)-lambda(1));
phi(:,length(lambda))=(p2(:,length(lambda))-p2(:,length(lambda)-1))/(lambda(length(lambda))-lambda(length(lambda)-1));
for j=2:length(lambda)-1
    phi(:,j)=(p2(:,j+1)-p2(:,j-1))/(lambda(j+1)-lambda(j-1));
    %plot(sigma,phi(:,j))
    %pause(.01)
end

psi(1,:)=(p2(2,:)-p2(1,:))/(sigma(2)-sigma(1));
psi(length(sigma),:)=(p2(length(sigma),:)-p2(length(sigma)-1,:))/(sigma(length(sigma))-sigma(length(sigma)-1));
for k=2:length(sigma)-1
    psi(k,:)=(p2(k+1,:)-p2(k-1,:))/(sigma(k+1)-sigma(k-1));
end

disp('Plotting lambda')
for j=1:length(lambda)
    plot(sigma,p2(:,j))
    pause(.1)
end
disp('finding Physical-1/5')
[LAM, Fgrid] = meshgrid(-lambda, m/(m+1)*sigma);
[~, intgrid] = meshgrid(-lambda, 1/2*m/(m+1)*(sigma.^2-sigma(1)^2));
disp('finding Physical-2/5')
u2   = psi./Fgrid;
eta2 = (phi-u2.^(2))/(2*g);
disp('finding Physical-3/5')
t2   = abs(((LAM-u2)/(alpha*g)));
x2   = (phi-u2.^(2)-intgrid)/(2*alpha*g);
disp('finding Physical-4/5')

eta2 = eta2(2:end-1,:);
t2   =   t2(2:end-1,:);
x2   =   x2(2:end-1,:);
u2   =   u2(2:end-1,:);


disp('finding Physical-5/5')
 
 
[x_lin t_lin eta_lin u_lin] = toConstantTime(x2,t2, 1:max(max(t2)) ,eta2, u2);

for i=1:length(t_lin(1,:))
    plot(x_lin(:,i),eta_lin(:,i));
    axis([min(min(x_lin)) max(max(x_lin)) min(min(eta_lin)) max(max(eta_lin)) ])
    pause(.1)
end

 for i=1:length(x_lin(1,:))
    plot(x_lin(:,i),eta_lin(:,i),'.r');
    hold on
    x=[-40 3];
    plot(x,alpha*x)
    plot(0,0,'^b')
    hold off
    axis([min(x) max(x) min(min(eta_lin)) max(max(eta_lin))])
    pause(.1)
end


end