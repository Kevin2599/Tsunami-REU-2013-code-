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
		else
			m.writer = VideoWriter(m.movieName);
			open(m.writer);
		end
	end

	figure('visible',readOption(options,'visible','on'));
	m.plotHandle = gcf();
	if isfield(options,'subplot')
		m.subplotDim = options.subplot;
		subplot(options.subplot(1), options.subplot(2), 1)
	end

	if m.saveMovie
		m.frameLoc = [m.loc '/frames-', num2str(m.plotHandle)];
		if inoctave()
			system(['rm -rf ' m.frameLoc]); 
			system(['mkdir -p ' m.frameLoc]);
		else
			% set(gcf,'Renderer','zbuffer');
		end
	end
end