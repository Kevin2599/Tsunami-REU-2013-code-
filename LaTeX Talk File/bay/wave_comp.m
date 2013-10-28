x=[1:-.01:-10];
xw=0:-.01:-10;
xw2=-.3:-.01:-10;
alpha=.01;
y=alpha*x;


subplot(2,1,1)
plot(x,y)
hold on
plot(x,0*x)
plot(xw,.005*sin(xw),'r')
xlabel('x')
ylabel('Wave Height')

subplot(2,1,2)

plot(x,y)
hold on
plot(x,0*x)
plot(xw2,.005*sin(xw2*8),'r')
xlabel('x')
ylabel('Wave Height')