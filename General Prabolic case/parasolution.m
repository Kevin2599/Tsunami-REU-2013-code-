function [Phi,x_lin t_lin eta_lin u_lin,besselpart]=parasolution()
    Num_Roots=30;
    m=2;

    lambda_S=0;
    lambda_E=100;
    %plotcoef=0;


    alpha=.05;
    g=9.81;
    maxx=.25*m/(m+1)/g/alpha;
    x=-[0:.0001*maxx:maxx/10 maxx/10:.01*maxx:9*maxx/10 9*maxx/10:.0001*maxx:maxx];

    sigma=sqrt(-x*alpha*g*(m+1)/m)*2;
    Sigma_L=max(sigma);


    lambda=(lambda_S:lambda_E/200:lambda_E);

    disp('bessel')
    Zeros_of_bessel=besselzero((1/m)+1,Num_Roots-1,1);
    Lambda_n=-([0 Zeros_of_bessel']./Sigma_L).^2;% build the eigenvalues from zeros of bessel
    %disp(Zeros_of_bessel)
    %{
    if(m==2)
        plot(tan(Sigma_L*sqrt(-Lambda_n))-Sigma_L*sqrt(-Lambda_n))%compares eigen values to expected ones for m=2, shows deviation from 0 for each eigen value.
        pause(1)
    end
    %}
    %find constants
    disp('constants')
    C_n=zeros(Num_Roots,1);
    D_n=zeros(Num_Roots,1);
    C = {'r','k','b',[.5 .6 .7],'m','r','c','g',[.8 .2 .6]};
    %bessel2=besselj(1/m,sqrt(-1*Lambda_n(2))*sigma)./sigma.^(1/m);
    for n=1:Num_Roots%need to make besselj better so that x can be 0 THIS MIGHT BE DONE
        
        besselpart=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
        if min(sigma)==0    
            besselpart(1)=(sqrt(-Lambda_n(n))/2)^(1/m)/gamma(1.+1/m);
        end
        
    %     normalized=besselpart./besselpart(1);
    %     
    %     if Lambda_n(n)==0
    %         normalized(:)=0;
    %     end
    %     disp('endif')
    %     %trapz(sigma,besselpart.*bessel2.*sigma.^(2/m+1))
    %     
    %     plot(sigma, normalized, 'color', C{n})
    %     legend(['1st (Eigenvalue= ' num2str(Lambda_n(1)) ')'], ['2nd (Eigenvalue= ' num2str(Lambda_n(2)) ')'] , ['3rd (Eigenvalue= ' num2str(Lambda_n(3)) ')'])
    %     xlabel('Sigma')
    %     ylabel('Eigenfunction Amplitude (Normalized)')
    %     disp('axis set')
    %     %axis([min(sigma) max(sigma) min(normalized)*1.2 1.2])
    %     pause(1)
    %     hold on
        
        C_n(n)=trapz(sigma,(sqrt(-1*Lambda_n(n))*2*g*eta_0(-sigma).*besselpart));
        D_n(n)=trapz(sigma(2:end),(besselpart(2:end).*cumsimps(sigma(2:end),m/(m+1).*sigma(2:end).*eta_0(-sigma(2:end)).*sqrt(g./(alpha*-x(2:end))))));
        %D_n(n)=0;
        
    end
    hold off





    %if plotcoef
    %    plot(C_n)
    %    hold on
    %    plot(D_n,'r')
    %    return;
    %end

    % NEED_ROOTS=max([abs(max(C_n(floor(Num_Roots*.9):end))),abs(min(C_n(floor(Num_Roots*.9):end))),abs(max(D_n(floor(Num_Roots*.9):end))),abs(min(D_n(floor(Num_Roots*.9):end)))]);%measure of how far from 0 last 10% of constants are
    % 
    % 
    % if NEED_ROOTS>10^-90
    %     disp('Need More Roots')
    %     plot(C_n)
    %     hold on
    %     plot(D_n,'r')
    %     hold off
    % 
    %     while NEED_ROOTS>10^-90 %adds more roots so that eigen expansion is more accurate if need be; make more restrictive for larger values of xmax.
    %         Num_Roots
    %         Num_Roots=Num_Roots+1;
    %         temp=besselzero((1/m)+1,Num_Roots-1,1);
    %         Lambda_n(Num_Roots)=-1*(temp(end)./Sigma_L).^2;%add the new lambda
    %         n=Num_Roots;
    %         besselpart=besselj(1/m,sqrt(-1*Lambda_n(n))*sigma)./sigma.^(1/m);
    %         C_n(n)=trapz(sigma,(sqrt(-1*Lambda_n(n))*2*g*eta_0(-sigma).*besselpart));
    %         D_n(n)=trapz(sigma,(besselpart.*cumsimps(sigma,m/(m+1).*sigma.*eta_0(-sigma).*sqrt(g./(alpha*-x)))));
    %         %D_n(n)=0;
    %         NEED_ROOTS=max(max(C_n(floor(Num_Roots*.9):end))-min(C_n(floor(Num_Roots*.9):end)),max(D_n(floor(Num_Roots*.9):end))-min(D_n(floor(Num_Roots*.9):end)));
    %         plot(C_n)%visual aid to see last portion of constants approaching 0 as new roots are added
    %         hold on
    %         plot(D_n,'r')
    %         hold off
    %         pause(1)
    %     end
    %     
    % if(m==2)%compares eigen values to expected ones for m=2, shows deviation from 0 for each eigen value
    %     plot(tan(Sigma_L*sqrt(-Lambda_n))-Sigma_L*sqrt(-Lambda_n))
    %     pause(1)
    % end
    % 
    % end



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


    disp('finding phi')


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
    disp('finding psi')
    psi(1,:)=(p2(2,:)-p2(1,:))/(sigma(2)-sigma(1));
    psi(length(sigma),:)=(p2(length(sigma),:)-p2(length(sigma)-1,:))/(sigma(length(sigma))-sigma(length(sigma)-1));
    for k=2:length(sigma)-1
            psi(k,:)=(p2(k+1,:)-p2(k-1,:))/(sigma(k+1)-sigma(k-1));
    end

    % for j=1:length(lambda)
    %     plot(sigma,p2(:,j))
    %     pause(.1)
    % end
    disp('final conversion')
    [LAM, Fgrid] = meshgrid(-lambda, m/(m+1)*sigma);
    [~, intgrid] = meshgrid(-lambda, 1/2*m/(m+1)*(sigma.^2-sigma(1)^2));

    u2   = psi./Fgrid;
    eta2 = (phi-u2.^(2))/(2*g);
    t2   = abs(((LAM-u2)/(alpha*g)));
    x2   = (phi-u2.^(2)-intgrid)/(2*alpha*g);

    eta2 = eta2(2:end-1,:);
    t2   =   t2(2:end-1,:);
    x2   =   x2(2:end-1,:);
    u2   =   u2(2:end-1,:);

    [x_lin t_lin eta_lin u_lin] = toConstantTime(x2,t2, 1:max(max(t2)) ,eta2, u2);

end