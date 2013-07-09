% function waterOutline = topViewOfWater(bayShape, waveX, waveHeight)
%
%  bay = {.left .right .height .slope}
%  waveX      (M x 1)
%  waveHeight (M x 1)
% ========
%  wateroutline (2M x 2) [x y]
% just do plot(wateroutline(:,1), wateroutline(:,2));

function [xs ys dys] = topViewOfWater(bay, waveX, waveHeight, magnification)
	waveHeight = waveHeight - bay.slope*waveX;

	outlineLeft  = interp1(bay.height, bay.left , waveHeight);
	outlineRight = interp1(bay.height, bay.right, waveHeight);

	origOutlineLeft  = interp1([bay.height; -1e20],[bay.left ;  bay.left(end)], -bay.slope*waveX);
	origOutlineRight = interp1([bay.height; -1e20],[bay.right; bay.right(end)], -bay.slope*waveX);
	
	dys = (outlineRight - outlineLeft) - (origOutlineRight - origOutlineLeft);
	dys = dys(end:-1:1);

	if exist('magnification') && magnification ~= 1
		outlineLeft  = (1-magnification) * origOutlineLeft  + magnification * outlineLeft;
		outlineRight = (1-magnification) * origOutlineRight + magnification * outlineRight;
	end

	xs = [waveX(end:-1:1); waveX];
	ys = [outlineLeft(end:-1:1); outlineRight];
end