% function waterOutline = topViewOfWater(bayShape, slope, waveX, waveHeight)
%
%  bayShape   (M x 3) [lY rY h]
%  slope      (1 x 1)
%  waveX      (M x 1)
%  waveHeight (M x 1)
% ========
%  wateroutline (2M x 2) [x y]
% just do plot(wateroutline(:,1), wateroutline(:,2));

function waterOutline = topViewOfWater(bayShape, slope, waveX, waveHeight)
	waveHeight = waveHeight - slope*waveX;

	bayLeft   = bayShape(:,1);
	bayRight  = bayShape(:,2);
	bayHeight = bayShape(:,3);

	outlineLeft  = interp1(bayHeight,bayLeft ,waveHeight);
	outlineRight = interp1(bayHeight,bayRight,waveHeight);

	waterOutline = [waveX(end:-1:1) outlineLeft(end:-1:1); waveX outlineRight];
end