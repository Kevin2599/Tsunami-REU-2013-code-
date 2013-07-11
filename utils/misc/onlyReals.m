function r = onlyReals(nums,cutoff=0.0001)
	r = [];
	for j=1:size(nums,2)
		is = find(abs(imag(nums(:,j))) < cutoff);
		for i=is'
			r(end+1,:) = [i, j, real(nums(i,j))];
		end
	end
end