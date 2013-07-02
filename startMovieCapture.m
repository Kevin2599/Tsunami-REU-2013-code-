function m = startMovieCapture(location)
	m.loc = location;
	m.frame = 0;
	if inoctave()
		system(sprintf('rm -rf %s; mkdir %s',location,location));
		clf; figure('visible','off');
	end
end