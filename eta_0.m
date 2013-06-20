function y=eta_0(x)
y=-9.0315e-4*exp(-1.5e-5*(1000+x).^2).*(1000+x);% gause wave
y(abs(y)<1e-3)=0;
y(1:1:length(y))=0;
end