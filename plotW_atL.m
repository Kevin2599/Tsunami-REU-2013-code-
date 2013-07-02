function [ll data] = zerosPowerSeries(f,max_precision)
	minmaxes = [0];
	curr_precision = 3;
	while curr_precision < max_precision
		[a b] = [f(0,curr_precision) f(0.01,curr_precision)];
end

function [x xx xxx] = diff2(f,x,dx=0.01)
	[a b c] = [f(x-dx) f(x) f(x+dx)];
	[x xx xxx] = [b (c-a)/(2*dx) (c - 2*b + a)/(dx^2)]