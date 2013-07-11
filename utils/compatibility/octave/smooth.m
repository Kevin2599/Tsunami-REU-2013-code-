% yy = smooth(y, span=5)
% smooth takes noisy data and smoothes it out
%
% Smooth is a Matlab function which has different smoothing methods;
%  this just implements the 'moving average' method.
% http://www.mathworks.com/help/curvefit/smooth.html

function yy = smooth(y, span=5)
	half = floor(span/2);
	if half*2 + 1 ~= span
		error('span must be odd: span=%d',span)
	end
	len = length(y);

	% filter is a rolling average, not a windowed average
	yy = filter(ones(1,span), span, y);

	% the middle bits of filter are a-okay
	yy(half+1:end-half) = yy(span:end);
	for i=1:half
		yy(i) = mean(y(1:i*2-1));
		yy(end-i+1) = mean(y(end-2*i+2:end));
	end

	%output must be a collumn vector
	if size(yy)(1) == 1
		yy = yy';
	end
end