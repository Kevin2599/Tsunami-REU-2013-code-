function z = zerosRange(f,noise_cutoff=.1,xs=0:.1:100)
	%initial pass to discover ranges to look for zeros
	ys = [];
	z = [];
	for x=xs;
		ys(end+1) = f(x);
	end

	% power series can get super noisy
	noise = filter(ones(1,10)/10,10, abs(diffn(ys,10)) );
	noise_over = find(noise > noise_cutoff);
	if length(noise_over) > 0
		xs = xs(1:noise_over(1));
		ys = ys(1:noise_over(1));
	end

	%use the search ranges to find the actual zeros
	idx = find(sign(ys(1:end-1).*ys(2:end)) < 0);
	for lr = [xs(idx) ; xs(idx+1)]
		z(end+1) = fzero(f,lr);
	end
end
function ys = diffn(ys,n)
	for i=1:n
		ys = diff(ys);
	end
end