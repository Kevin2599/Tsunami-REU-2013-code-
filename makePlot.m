function m = makePlot(varargin)
	m = readOptions(varargin);

	m.saveMovie = readOption(m,'saveMovie',false);

	if m.saveMovie
		m.movieName = readOption(m,'movieName','octave.avi');
		m.loc = readOption(m,'movieLocation','octaveMovies');
		m.visible = readOption(m,'visible','off')
		m.closeOnFinish = true;
		m.frame = 0;
		if inoctave()
			graphics_toolkit gnuplot;
		else
			m.writer = VideoWriter(m.movieName);
			open(m.writer);
		end
	end

	figure('visible',readOption(m,'visible','on'));
	m.plotHandle = gcf();
	if isfield(m,'subplot')
		m.subplotDim = m.subplot;
		subplot(m.subplot(1), m.subplot(2), 1)
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