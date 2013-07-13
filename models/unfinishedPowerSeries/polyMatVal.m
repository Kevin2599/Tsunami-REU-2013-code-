function y = polyMatVal(C,x)
	for i=1:size(C,2)
		y(:,i) = polyval(C(:,i),x);
	end
end