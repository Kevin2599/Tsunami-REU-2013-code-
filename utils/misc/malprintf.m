function malprintf(varargin)
	fprintf(varargin{:});
	flushConsole();
end