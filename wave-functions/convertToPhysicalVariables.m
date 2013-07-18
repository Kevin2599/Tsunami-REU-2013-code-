function vars = convertToPhysicalVariables(Phi, Psi, lambda, F, intF, options)
	g = options.g;
	alpha = options.bath.slope;


	if size(Phi,2) ~= length(lambda) % basic error checking
		Phi    =    Phi'; % sigma x lambda
		Psi    =    Psi'; % sigma x lambda
	end

	% Data Needed to convert both exact and aprox data
	[LAM, Fgrid] = meshgrid(-lambda, F);
	[~, intgrid] = meshgrid(-lambda, intF);

	% Convert Aprox.
	vars.u   = Psi./Fgrid;
	vars.t   = abs((LAM-vars.u)/(alpha*g));
	vars.eta = (Phi-vars.u.^(2))/(2*g);
	vars.x   = (Phi-vars.u.^(2)-intgrid)/(2*alpha*g);

	%[J, UL, US]=Jacobian(F,g,alpha,u,sigma,lambda,dsigma,dlambda);

	% Test the data to see if it breaks
	[~,I]=sort(vars.t,2);
	brokeat=0;
	for j=1:length(vars.t(1,:))
		if I(2,j)~=j
			brokeat=j;
			println(['Broke at ' num2str(brokeat)])
			if options.trimAtBreak
				vars = applyFunToStruct(@(v) v(:, 1:j-1) , vars);
			end
			break
		end
	end
end