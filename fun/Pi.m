n=10000000;
in=0;
for i=1:n
    x=2*rand(2,1)-1;
    
    if (x(1)^2+x(2)^2<=1)
        in=in+1;
    end
    
end
4*in/n-pi
