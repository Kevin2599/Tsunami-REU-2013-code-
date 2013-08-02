function flushConsole()
	if inoctave()
		fflush(stdout);
	else
		drawnow('update');
	end
end