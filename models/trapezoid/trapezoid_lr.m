function [l r] = trapezoid_lr(height, width, slope)
	r = (height >= 0) .* (height.*slope + width/2);
	l = -r;
end