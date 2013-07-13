function y = initialHump(amplitude, offset, width, x)
	y = amplitude * exp( - ((x - offset)/width).^2);
end