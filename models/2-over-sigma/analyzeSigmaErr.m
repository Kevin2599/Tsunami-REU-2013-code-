runup = zeros(length(results),4);

function [m idx] = matMax(mat)
	[mx i] = max(mat);
	[m  j] = max(mx);
	idx = {i(j) j};
end

for i=1:length(results)
	r = results{i};
	[runup(i,1) runup(i,3)] = min(r.eta(1,:));
	[runup(i,2) runup(i,4)] = max(r.eta(1,:));
	runup(i,3) = r.t(1,runup(i,3));
	runup(i,4) = r.t(1,runup(i,4));
end

semilogy(abs(runup));