function [J, U_lambda, U_sigma] = Jacobian(F,g,alpha,U,sigma,lambda,dsigma,dlambda)

for i=1:size(U,1)
    U_lambda(i,1)=(U(i,1)-U(i,2))/dlambda;
    for j=2:size(U,2)-1
        U_lambda(i,j)=(U(i,j+1)-U(i,j-1))/2/dlambda;
    end
    U_lambda(i,size(U,2))=(U(i,end)-U(i,end-1))/dlambda;
end
%U_lambda(sigma,lambda)

for j=1:size(U,2)
    U_sigma(1,j)=(U(1,j)-U(2,j))/dsigma;
    for i=2:size(U,1)-1
        U_sigma(i,j)=(U(i+1,j)-U(i-1,j))/2/dsigma;
    end
    U_sigma(size(U,1),j)=(U(end,j)-U(end-1,j))/dsigma;
end
%U_sigma(sigma,lambda)

for i=1:length(lambda)
    for j=1:length(sigma)
        %J(j,i)=F(j)/(4*g^2*alpha^2)*((1-U_lambda(j,i))^2-U_sigma(j,i)^2);
        J(j,i)=(1-U_lambda(j,i))^2-U_sigma(j,i)^2;
        if (J(j,i)==0)
            disp('J=0!!!')
        end
    end
end
end