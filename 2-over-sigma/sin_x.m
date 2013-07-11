function y = sin_x(x,n=5)
	y = ones(size(x));
	for i=1:n
		y = y + (-1)^i * x.^(2*i) ./ factorial(2*i + 1);
	end
end