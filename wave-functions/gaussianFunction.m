function y = gaussianFunction(offset, width, x, derivative)
	y = exp( - ((x - offset)./width).^2);

	if exist('derivative')
		switch derivative
		case 0
			y=y;
	 	case 1
			y = y.* -2 .*((x - offset)./width);
		end
	end
end