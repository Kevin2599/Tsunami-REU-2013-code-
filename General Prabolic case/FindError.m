error=NaN*ones(length(sigma)/10, 1);
for i=1:10:length(sigma)
     p=0;
    for n=1:length(C_n)
        p=p+C_n(n)*(sqrt(-1*Lambda_n(n)))*besselpart(n,i);
    
    end
        error(1+(i-1)/10)=p-phi_0(sigma(i));
end