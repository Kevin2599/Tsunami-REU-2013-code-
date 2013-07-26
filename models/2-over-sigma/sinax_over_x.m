function y = sinax_over_x(a,x,derivative)
	if(length(a) == 1)
		a = a * ones(size(x));
	elseif(length(x) == 1)
		x = x * ones(size(a));
	end

	if ~exist('derivative')
		derivative = 0;
	end

	ax = a.*x;

	% with cutoff = .1, |f' - f|/|1-f| ~ 1.2e-7
	close0 = ax <.1;
	far0   = ~close0;

	y = zeros(size(ax));

	switch derivative
		case 0
			y(close0) = a(close0)    .* (1 - ax(close0).^2/6 + (ax(close0).^4)/120);
			y(far0  ) = sin(ax(far0)) ./ x(far0);
		case 1
			y(close0) = a(close0).^2 .* (  - ax(close0)   /3 + (ax(close0).^3)/30);
			y(far0  ) = (ax(far0) .* cos(ax(far0)) - sin(ax(far0))) ./ x(far0).^2;
	end

	% d/dx(f) ~ a^2 [-ax/3 + (ax)^3 / 30 - (ax)^5 / 840]
	% df < | (ax)^5 / 840 | / | -ax/3 + (ax)^3 / 30 |
	%    x<-2 (7 df+sqrt(7) sqrt(df (7 df+10)))
	%    with cutoff = .1, |f' - f|/f ~ 3.5e-7
end