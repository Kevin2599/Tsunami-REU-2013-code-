function m = makePlot(plotID,varargin)
	m = readOptions(varargin);

	global numPlots = 0;
	numPlots = numPlots+1;

	m.plotID = plotID;
	m.plotHandle = -1;
	m.closeOnFinish = false;
	m.frame = 0;
	m.writer = [];
	m.frameLoc = [];

	m.saveMovie = readOption(m,'saveMovie',false);

	if m.saveMovie
		m.movieName     = readOption(m,'movieName','octave.avi');
		m.loc           = readOption(m,'movieLocation','octaveMovies');
		m.visible       = readOption(m,'visible','off');
		m.closeOnFinish = true;
		if inoctave()
			graphics_toolkit gnuplot;
		else
			m.writer = VideoWriter(m.movieName);
			open(m.writer);
		end

		m.frameLoc = [m.loc '/frames-', num2str(numPlots)];
		if inoctave()
			system(['rm -rf ' m.frameLoc]); 
			system(['mkdir -p ' m.frameLoc]);
		else
			% set(gcf,'Renderer','zbuffer');
		end
	end
end