function printf(varargin)
	fprintf(varargin{:});
	if inoctave()
		fflush(stdout);
	else
		drawnow('update');
	end
end