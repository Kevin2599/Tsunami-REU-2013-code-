function [ Besselsig BesselsigNorm eta x u t] = BesselGenerator( n,m, alpha )%Root, power of bay, slope of bay
close all

tic
Num_Roots=n;


lambda_S=0;
lambda_E=50;

plotb=0;
%put numbers for our problem

Sigma_L=80;
g=1;%9.81;
%sigma=[0:.001:10 10.001:.01:Sigma_L];
sigma=0:.1:Sigma_L;
% sigma=zeros(200,1);
% sigma(1)=Sigma_L*1.1;
% dx=Sigma_L/10;
% factor=1.01;
% for i=2:200
%     dx=dx*factor;
%     sigma(i)=sigma(i-1)+dx;
% end
% sigma=sigma'*Sigma_L/max(sigma);
% plot(sigma,'.')


lambda=(lambda_S:.1:lambda_E);
disp('Zeros of bessel')
Zeros_of_bessel=besselzero((1/m)+1,Num_Roots-1,1);
Lambda_n=-( Zeros_of_bessel(end)./Sigma_L).^2;% build the eigenvalues from zeros of bessel


%find constants
disp('Finding Consts using Eigenfunctions')



%From the inner product
% besselpart=zeros(Num_Roots,length(sigma));
% for n=1:Num_Roots%need to make besselj better so that x can be 0
%
%     besselpart(n,:)=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
%     %fix the zero point
%     besselpart(n,1)=besseljo(1/m,Lambda_n(n),0);
%
%
%     [phi pot,A, p, sigma_0]=phi_0(sigma);
%     C_n(n)=trapz(sigma,(sqrt(-1*Lambda_n(n))*phi.*besselpart(n,:)));%(-Lambda_n(n))^.25/sqrt(2/pi).*
% %     f=@(xt) (sqrt(-1*Lambda_n(n))).*besselj(1/m,sqrt(-1*Lambda_n(n))*xt)./xt.^(1/m).*...
% %         (-4*A*xt.^(-1).*((xt-sigma_0)/p^2.*exp(-1*((xt-sigma_0)/p).^2)+(xt+sigma_0)/p^2.*exp(-((xt+sigma_0)/p).^2)))';
% %     C_n(n)=integral(f,0,Sigma_L);
%     D_n(n)=pot;
%     if (plotb)
%         if(n==Num_Roots)
%             figure(1)
%             plot(sigma,besselpart(n,:),'.b')
%             hold on
%             plot(sigma(1:end-1),diff(besselpart(n,:))./diff(sigma),'.r','MarkerSize',1)
%             plot(sigma,0.*sigma,'k')
%             legend('Eig(sigma)','Eigp(sigma)');
%             hold off
%
% %             figure(2)
% %             plot(x,besselpart(n,:),'.b')
% %             hold on
% %             plot(x(1:end-1),diff(besselpart(n,:))./diff(x),'.r','MarkerSize',1)
% %             plot(x,0.*x,'k')
% %             legend('Eig(x)','Eigp(x)');
% %             set(gca,'xdir','reverse')
% %             hold off
%         end
%     end
% end


% using the matrix.
besselpart(1,:)=besselj(1/m,sqrt(-1*Lambda_n)*sigma)./sigma.^(1/m);
%fix the zero point
besselpart(1,1)=besseljo(1/m,Lambda_n,0);
 Besselsig=besselpart;
 BesselsigNorm=Besselsig/max(Besselsig);
 Amatrix(:,1)=(sqrt(-1*Lambda_n))*besselpart(1,:);
%[phi pot,A, p, sigma_0]=phi_0(sigma);
pot=0;
phi=besselpart(1,:);
C_n=[ phi/Amatrix'];
D_n=C_n*pot;
toc


figure(1)
plot(sigma,Besselsig)
xlabel('sigma')
ylabel('heigth')
title('Besselsig')
figure(2)
plot(sigma,BesselsigNorm)
xlabel('sigma')
ylabel('heigth')
title('BesselsigNorm')
disp('Building function')
p2= zeros(length(sigma),length(lambda));

for j=1:length(lambda)
    p2(:,j)=besselpart(1,:).*(C_n(1).*sin(sqrt(-1*Lambda_n(1))*lambda(j))+D_n(1).*cos(sqrt(-1*Lambda_n(1))*lambda(j)));
end


disp('finding phi and psi...')


phi=zeros(length(sigma),length(lambda)); %had length(lambda)-1 before, not necessary

phi(:,1)=(p2(:,2)-p2(:,1))/(lambda(2)-lambda(1));
phi(:,length(lambda))=(p2(:,length(lambda))-p2(:,length(lambda)-1))/(lambda(length(lambda))-lambda(length(lambda)-1));
for j=2:length(lambda)-1
    phi(:,j)=(p2(:,j+1)-p2(:,j-1))/(lambda(j+1)-lambda(j-1));
end


psi=zeros(length(sigma),length(lambda));

psi(1,:)=(p2(2,:)-p2(1,:))/(sigma(2)-sigma(1));
psi(length(sigma),:)=(p2(length(sigma),:)-p2(length(sigma)-1,:))/(sigma(length(sigma))-sigma(length(sigma)-1));
for k=2:length(sigma)-1
    psi(k,:)=(p2(k+1,:)-p2(k-1,:))/(sigma(k+1)-sigma(k-1));
end


phi=m/(m+1)*phi;
psi=m/(m+1)*psi;

disp('final conversion')
[LAM, Fgrid] = meshgrid(-lambda, m/(m+1)*sigma);
[~, intgrid] = meshgrid(-lambda, m/2*(1/(m+1))*(sigma.^2));



u   = psi./Fgrid;
eta = (phi-u.^(2))/(2*g);
t   = abs(((LAM-u)/(alpha*g)));
x   = (phi-u.^(2)-intgrid)/(2*alpha*g);

eta = eta(2:end-1,:); %WHY CUT OUT ENDPOINTS?
t   =   t(2:end-1,:);
x  =   x(2:end-1,:);
u   =   u(2:end-1,:);



%[x_lin t_lin eta_lin u_lin] = toConstantTime(x2,t2, 1:max(max(t2)) ,eta2, u2);
%WHY IS THIS COMMENTED?




%Graph%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc



%final plot our model

figure(3)
for i=2:10:length(p2(1,:))
    plot(sigma,p2(:,i))
    hold on
    plot([max(max(sigma)) min(min(sigma))] ,[0 0],'k')
    legend('error','zero')
    xlabel('Sigma')
    pause(.01)
    hold off
end


toc
figure(4)
scale=50;
for i=1:10:length(lambda)
    plot(x(:,i),eta(:,i),'.b')
    hold on
    %plot(x2(:,i),u2(:,i),'.r')
    plot([max(max(x)) -scale*max(max(x))],alpha*[max(max(x)) -scale*max(max(x))],'k') %notation?
    plot(0,0,'^k')
    set(gca,'xdir','reverse')
    xlabel('x')
    legend('eta(x)','error','bay')
    axis([-scale*max(max(x)) max(max(x)) scale/5*min(min(eta)) scale/5*max(max(eta))])
    title(['t = ' num2str(t(2,i))]);
    hold off
    pause(.1)
end





% aviobj = VideoWriter('Runup.avi');
% fig=figure(4);
% scale=25;
% open(aviobj)
% for i=100:1:193
%     plot(x2(:,i),eta2(:,i),'.b')
%     hold on
%      plot(x2(:,i),abs(eta2(:,i)-eta2E(:,i)),'r')
%     %plot(x2(:,i),u2(:,i),'.r')
%     plot([max(max(x2)) -scale*max(max(x2))],alpha*[max(max(x2)) -scale*max(max(x2))],'k') %notation?
%     plot(0,0,'^k')
%     set(gca,'xdir','reverse')
%     legend('eta(x)','error','bay')
%     axis([-scale*max(max(x2)) max(max(x2)) scale/5*min(min(eta2)) scale/5*max(max(eta2))])
%     title(['t = ' num2str(t2(2,i))]);
%     xlabel('x')
%     ylabel('eta')
%     hold off
%
%     pause(.1)
%      F = getframe(fig);
%     writeVideo(aviobj,F);
% end
% close(fig);
%  close(aviobj);
%
%
%
%












end
