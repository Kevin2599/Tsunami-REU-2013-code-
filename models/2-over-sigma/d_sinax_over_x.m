function y = d_sinax_over_x(a,x)

	if(length(a) == 1)
		a = a * ones(size(x));
	elseif(length(x) == 1)
		x = x * ones(size(a));
	end

	ax = a.*x;

	close0 = ax < .1;
	far0   = ~close0;

	y = zeros(size(ax));

end