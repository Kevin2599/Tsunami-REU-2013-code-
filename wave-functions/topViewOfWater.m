% function [xs ys dyw] = topViewOfWater(bayShape, waveX, waveHeight)
%
%  bay = {.lr .slope}
%  waveX      (M x 1)
%  waveHeight (M x 1)
% ========
% just do plot(xs,ys);

function [xs ys dys] = topViewOfWater(bay, waveX, waveHeight, magnification)
	initialHeight =  - bay.slope*waveX;

	[outlineLeft outlineRight] = bay.lr(waveHeight + initialHeight, bay);

	if exist('magnification') && magnification ~= 1
		[origOutlineLeft origOutlineRight] = bay.lr(initialHeight, bay);
		
		dys = (outlineRight - outlineLeft) - (origOutlineRight - origOutlineLeft);
		dys = dys(end:-1:1);

		outlineLeft  = (1-magnification) * origOutlineLeft  + magnification * outlineLeft;
		outlineRight = (1-magnification) * origOutlineRight + magnification * outlineRight;
	else
		dys = [];
	end

	xs = [waveX(end:-1:1); waveX];
	ys = [outlineLeft(end:-1:1); outlineRight];
end