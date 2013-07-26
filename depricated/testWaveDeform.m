figure(1); clf; hold on
x_axis = [-100 10];
y_axis = [-1 1];
mag = 5;
alpha = .05;


function y = logistic(x,x0,a)
    y = 1 ./ (1 + exp((x0-x).*a));
end

function y = magnifyLogistic(ys,xs,mag,alpha)
	y = ys .* (1 + mag .* logistic(-xs,10,1));
end

function y = magnifyBlend(ys, xs, mag, alpha)
	y1 = mag*ys; y2 = ys;
	x1 = -abs(2*mag*ys/alpha) - 10; x2 = ys./alpha;
	m = (y2 - y1)./(x2 - x1);

	left = mag*ys; right = m .* (xs-x1) + y1;

	a = logistic(xs - x1./2,-10,abs(xs-x1)./50);

	y = (1-a) .*left + a .*right;
end

function y = magnify(ys, xs, mag, alpha)
	y1 = mag*ys; y2 = ys;
	x1 = -abs(2*mag*ys/alpha) - .1e-17; x2 = ys./alpha;
	m = (y2 - y1)./(x2 - x1);

	plot(x1,y1,'^b');
	plot(x2,y2,'^b');

	% this is the hyperbolic y=1/x that has been skewed by [1 -1/m; 0 1]
	b = -m.*(xs - alpha.*ys + ((ys.*(1-mag)).^2 - abs(m))./(m .* ys .* (1-mag)));
	c = -abs(m);
	% c = -ys.^2 .*(1 - mag./alpha);

	y = (-b - sign(ys).* sqrt(b.^2 - 4.*c))./2 + y1;
end

xs = linspace(x_axis(1),x_axis(2),50);

for mag = (1.01:.02:1.5).^4

	hold off; plot(x_axis,alpha*x_axis*.5,'k');
	hold on;  plot(x_axis,alpha*x_axis*-.5,'k');

	plot(x_axis,alpha*x_axis);
	plot(0,0,'^b');
	for eta = linspace(y_axis(1)/mag, y_axis(2)/mag, 17)

		% plot(xs,magnifyBlend(eta.*ones(size(xs)), xs, mag, alpha),'g');
		plot(xs,magnify(eta.*ones(size(xs)), xs, mag, alpha),'r');
	end
	xlim(x_axis);
	ylim(y_axis);
	title(['Magnification ' num2str(mag)])
	drawnow();
end

