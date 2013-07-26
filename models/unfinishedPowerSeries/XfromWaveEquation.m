function A = XfromWaveEquation(W,n0,L,zeroCutoff = 1e-3)
	X = @(lambda,precision) polyval (solve_Xxx_WXx_X (W(1:precision), n0, lambda),L);

	zBig = zerosRange(@(x) X(x,length(W)));
	zSmall = zerosRange(@(x) X(x,length(W) -1 + n0));

	err = abs((zBig -zSmall)./zBig)
	lambda = zBig(find( err < zeroCutoff ))
	A=[];
	for l=lambda
		A(:,end+1) = solve_Xxx_WXx_X (W, n0, l);
	end
end