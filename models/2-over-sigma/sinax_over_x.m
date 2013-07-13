function y = sinax_over_x(a,x)
	if(length(a) == 1)
		a = a * ones(size(x));
	elseif(length(x) == 1)
		x = x * ones(size(a));
	end

	ax = a.*x;

	% with cutoff = .1, |f' - f|/|1-f| ~ 1.2e-7
	close0 = ax <.1;
	far0   = ~close0;

	y = zeros(size(ax));

	y(close0) = a(close0) .* (1 - ax(close0).^2/6 + (ax(close0).^4)/30);
	y(far0  ) = sin(ax(far0)) ./ x(far0);
end