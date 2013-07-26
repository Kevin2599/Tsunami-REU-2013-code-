clear all
tic
Num_Roots=90;
m=2;

lambda_S=0;
lambda_E=50;

plotb=0;
%put numbers for our problem
alpha=.05;
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
Lambda_n=-([0 Zeros_of_bessel']./Sigma_L).^2;% build the eigenvalues from zeros of bessel
if(m==2&&randi(100)==100)%plots about 1% of the time lul
    plot(tan(Sigma_L*sqrt(-Lambda_n))-Sigma_L*sqrt(-Lambda_n))%compares eigen values to expected ones for m=2, shows deviation from 0 for each eigen value.
    pause(1)
end

%find constants
disp('Finding Consts using Eigenfunctions')
C_n=zeros(Num_Roots,1);
D_n=zeros(Num_Roots,1);



%From the inner product
besselpart=zeros(Num_Roots,length(sigma));
for n=1:Num_Roots%need to make besselj better so that x can be 0

    besselpart(n,:)=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
    %fix the zero point
    besselpart(n,1)=besseljo(1/m,Lambda_n(n),0);


    [phi pot,A, p, sigma_0]=phi_0(sigma);
    C_n(n)=trapz(sigma,(sqrt(-1*Lambda_n(n))*phi.*besselpart(n,:)));%(-Lambda_n(n))^.25/sqrt(2/pi).*
%     f=@(xt) (sqrt(-1*Lambda_n(n))).*besselj(1/m,sqrt(-1*Lambda_n(n))*xt)./xt.^(1/m).*...
%         (-4*A*xt.^(-1).*((xt-sigma_0)/p^2.*exp(-1*((xt-sigma_0)/p).^2)+(xt+sigma_0)/p^2.*exp(-((xt+sigma_0)/p).^2)))';
%     C_n(n)=integral(f,0,Sigma_L);
    D_n(n)=pot;
    if (plotb)
        if(n==Num_Roots)
            figure(1)
            plot(sigma,besselpart(n,:),'.b')
            hold on
            plot(sigma(1:end-1),diff(besselpart(n,:))./diff(sigma),'.r','MarkerSize',1)
            plot(sigma,0.*sigma,'k')
            legend('Eig(sigma)','Eigp(sigma)');
            hold off

%             figure(2)
%             plot(x,besselpart(n,:),'.b')
%             hold on
%             plot(x(1:end-1),diff(besselpart(n,:))./diff(x),'.r','MarkerSize',1)
%             plot(x,0.*x,'k')
%             legend('Eig(x)','Eigp(x)');
%             set(gca,'xdir','reverse')
%             hold off
        end
    end
end


% using the matrix.
% Amatrix=zeros(length(sigma),Num_Roots-1);
% besselpart=zeros(Num_Roots,length(sigma));
% for n=2:Num_Roots
%     besselpart(n,:)=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
%     %fix the zero point
%     besselpart(n,1)=besseljo(1/m,Lambda_n(n),0);
%     
%     Amatrix(:,n-1)=(sqrt(-1*Lambda_n(n)))*besselpart(n,:);
% end
% [phi pot,A, p, sigma_0]=phi_0(sigma);
% C_n=[0 phi/Amatrix'];
% D_n=C_n*pot;


figure(3)
plot(C_n,'.-b')
hold on
plot(D_n,'.-r')
plot(1:length(D_n),0.*D_n,'k')
legend('C_n','D_n');
xlabel('n')
ylabel('constants')
hold off
pause(1)

disp('Checking convergence of constants')
if abs(max(C_n))+abs(min(C_n))+abs(max(D_n))+abs(min(D_n))==0
    disp('All constants are zero')
    return;
end

NEED_ROOTS=max([abs(max(C_n(floor(Num_Roots*.9):end))),abs(min(C_n(floor(Num_Roots*.9):end))),abs(max(D_n(floor(Num_Roots*.9):end))),abs(min(D_n(floor(Num_Roots*.9):end)))]);%measure of how far from 0 last 10% of constants are
if NEED_ROOTS>.01
    disp('May Need more roots for accuracy')
end



disp('Building function from basis')
Phi= zeros(length(sigma),length(lambda));

for j=1:length(lambda)
    Phi(:,j)=besselpart(1,:).*(C_n(1).*sin(sqrt(-1*Lambda_n(1))*lambda(j))+D_n(1).*cos(sqrt(-1*Lambda_n(1))*lambda(j)));
end

for n=2:Num_Roots
    for j=1:length(lambda)
        Phi(:,j)=Phi(:,j)+(besselpart(n,:).*(C_n(n).*sin(sqrt(-1*Lambda_n(n))*lambda(j))+D_n(n).*cos(sqrt(-1*Lambda_n(n))*lambda(j))))';
    end
    if( mod(floor((Num_Roots-n)/Num_Roots)*100,2)==0) %why doesn't this work correctly?
        disp([num2str(100-floor((Num_Roots-n)/Num_Roots*100)) '%'])
    end
end
p2(:,:)=Phi(:,:);

clear Phi
close all
if plotb
    figure(4)
    for i=1:length(p2(1,:))
        plot(sigma,p2(:,i))
        legend('Potential')
        pause(.01)
    end
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

if plotb
    for K=1:length(lambda)
        plot(sigma,psi(:,K),'r')
        hold on
        plot(sigma,phi(:,K))
        legend('psi','phi');
        hold off
        pause(.01)
    end
end
if m==2
    %lambda=-lambda;
    [LAM, SIG] = meshgrid(-lambda, sigma);
    
    % Find exact phi and psi.
    phiE = A./ SIG .* ( -2.*(SIG+LAM-sigma_0)/(p^2).* exp(-(((SIG+LAM-sigma_0)/p).^2)) -2.*(SIG-LAM-sigma_0)/(p^2).* exp(-(((SIG-LAM-sigma_0)/p).^2)) -2.*(SIG+LAM+sigma_0)/(p^2).* exp(-(((SIG+LAM+sigma_0)/p).^2))-2.*(SIG-LAM+sigma_0)/(p^2).* exp(-(((SIG-LAM+sigma_0)/p).^2)));
    psiE = -A./ SIG .* ( -(1./SIG + 2*(SIG + LAM - sigma_0)/(p^2)).* exp(-(((SIG + LAM - sigma_0)/p).^2)) + (1./SIG + 2*(SIG - LAM - sigma_0)/(p^2)).* exp(-(((SIG - LAM - sigma_0)/p).^2)) - (1./SIG + 2*(SIG + LAM + sigma_0)/(p^2)).* exp(-(((SIG + LAM + sigma_0)/p).^2)) + (1./SIG + 2*(SIG - LAM + sigma_0)/(p^2)).* exp(-(((SIG - LAM + sigma_0)/p).^2)));
    phiE=m/(m+1)*phiE;
    psiE=m/(m+1)*psiE;
end
phi=m/(m+1)*phi;
psi=m/(m+1)*psi;

disp('final conversion')
[LAM, Fgrid] = meshgrid(-lambda, m/(m+1)*sigma);
[~, intgrid] = meshgrid(-lambda, m/2*(1/(m+1))*(sigma.^2));



u2   = psi./Fgrid;
eta2 = (phi-u2.^(2))/(2*g);
t2   = abs(((LAM-u2)/(alpha*g)));
x2   = (phi-u2.^(2)-intgrid)/(2*alpha*g);

eta2 = eta2(2:end-1,:); %WHY CUT OUT ENDPOINTS?
t2   =   t2(2:end-1,:);
x2   =   x2(2:end-1,:);
u2   =   u2(2:end-1,:);


if m==2
    
    u2E   = psiE./Fgrid;
    eta2E = (phiE-u2E.^(2))/(2*g);
    t2E   = abs(((LAM-u2E)/(alpha*g)));
    x2E   = (phiE-u2E.^(2)-intgrid)/(2*alpha*g);
    
    eta2E = eta2E(2:end-1,:); %WHY CUT OUT ENDPOINTS?
    t2E   =   t2E(2:end-1,:);
    x2E   =   x2E(2:end-1,:);
    u2E   =   u2E(2:end-1,:);
end
%[x_lin t_lin eta_lin u_lin] = toConstantTime(x2,t2, 1:max(max(t2)) ,eta2, u2);
%WHY IS THIS COMMENTED?




%Graph%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

if m==2
    PPhi=zeros(length(sigma),length(lambda));
    for n=1:length(sigma)
        for j=1:length(lambda)
            PPhi(n,j)=A/sigma(n)*[exp(-((sigma(n)+lambda(j)-sigma_0)/p)^2)-exp(-((sigma(n)-lambda(j)-sigma_0)/p)^2)+exp(-((sigma(n)+lambda(j)+sigma_0)/p)^2)-exp(-((sigma(n)-lambda(j)+sigma_0)/p)^2)];
        end
    end
    
    
    figure(10)
    for i=2:length(p2(1,:))
        plot(sigma,abs(p2(:,i)-PPhi(:,i)))
        hold on
        plot(sigma,PPhi(:,i),'.r')
        plot([max(max(sigma)) min(min(sigma))] ,[0 0],'k')
        legend('error','pelinovsky','zero')
        pause(.01)
        hold off
    end
    
    [vm,im]=min(min(PPhi));
    [v,i]=max(max(PPhi));
    maxRunErr=abs(p2(2,i)-PPhi(2,i))/PPhi(2,i)
    minRunErr=abs(p2(2,im)-PPhi(2,im))/PPhi(2,im)
    
    % valud with few non zero points
    err=abs(p2(2,:)-PPhi(2,:))./PPhi(2,:);
    figure(5)
    loglog(x2(2,:),err)
end

pause(3)
%final plot our model
if m~=2
    figure(10)
    for i=2:length(p2(1,:))
        plot(sigma,p2(:,i))
        hold on
        plot([max(max(sigma)) min(min(sigma))] ,[0 0],'k')
        legend('error','zero')
        pause(.01)
        hold off
    end
end
if m==2
    toc
    figure(4)
    scale=25;
    for i=100:1:193
        plot(x2(:,i),eta2(:,i),'.b')
        hold on
        plot(x2(:,i),abs(eta2(:,i)-eta2E(:,i)),'r')
        %plot(x2(:,i),u2(:,i),'.r')
        plot([max(max(x2)) -scale*max(max(x2))],alpha*[max(max(x2)) -scale*max(max(x2))],'k') %notation?
        plot(0,0,'^k')
        set(gca,'xdir','reverse')
        legend('eta(x)','error','bay')
        axis([-scale*max(max(x2)) max(max(x2)) scale/5*min(min(eta2)) scale/5*max(max(eta2))])
        title(['t = ' num2str(t2(2,i))]);
        hold off
        pause(.1)
    end
end
if m~=2
    toc
    figure(4)
    scale=25;
    for i=1:1:length(lambda)
        plot(x2(:,i),eta2(:,i),'.b')
        hold on
        %plot(x2(:,i),u2(:,i),'.r')
        plot([max(max(x2)) -scale*max(max(x2))],alpha*[max(max(x2)) -scale*max(max(x2))],'k') %notation?
        plot(0,0,'^k')
        set(gca,'xdir','reverse')
        legend('eta(x)','error','bay')
        axis([-scale*max(max(x2)) max(max(x2)) scale/5*min(min(eta2)) scale/5*max(max(eta2))])
        title(['t = ' num2str(t2(2,i))]);
        hold off
        pause(.1)
    end
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










