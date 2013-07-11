function A = diff2Mat(m)
	A = (diag(ones(1,m-1),1) - 2*diag(ones(1,m)) + diag(ones(1,m-1),-1));
	A(1  ,:) = 2*A(2    ,:) - A(3    ,:);
	A(end,:) = 2*A(end-1,:) - A(end-2,:);
end