function m = makePlot(varargin)
	options = readOptions(varargin);

	m.saveMovie = readOption(options,'saveMovie',false);

	if m.saveMovie
		m.movieName = readOption(options,'movieName','octave.avi');
		m.loc = readOption(options,'movieLocation','octaveMovies');
		m.visible = readOption(options,'visible','off')
		m.closeOnFinish = true;
		m.frame = 0;
		if inoctave()
			graphics_toolkit gnuplot;
		end
	end

	figure('visible',readOption(options,'visible','on'));
	m.plotHandle = gcf();

	if m.saveMovie
		m.frameLoc = [m.loc '/frames-', num2str(m.plotHandle)];
		if inoctave()
			system(['rm -rfp ' m.frameLoc]); 
			system(['mkdir ' m.frameLoc]);
		end
	end
end