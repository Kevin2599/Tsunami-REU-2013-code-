%

function X = solve_Xxx_WXx_X(W,n0,lambda)
	if n0 > 0
		error('n0 can''t be 1 or greater');
	end
	W = reshape(W, 1, []);
	W = W(end:-1:1);	

	n = length(W)+1;
	M = zeros(n);
	x0 = 2 - n0;

	M(1,1) = 1;

	if n0 == 0
		M = M + diag(0:n-2+n0,1);
	else
		M = M + diag([0 0:n-1+n0],n0+1);
	end
	for i=2:n
		M(i,2:i) = M(i,2:i) + W(end-(i-2):end);
	end

	X = zeros(n,length(lambda));
	for i=1:length(lambda)
		X(:,i) = (M + diag(lambda(i)./[1 1:n-x0] , 1-x0)) \ eye(n,1);
		X(2:end,i) = X(2:end,i) ./ (1:n-1)';
	end
	X = X(end:-1:1,:);
end