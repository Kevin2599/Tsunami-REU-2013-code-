function m = drawPlot(m)
	if m.saveMovie
		if inoctave()
			fileName = sprintf('%s/%05d.png',m.frameLoc,m.frame);
			print(fileName);
			m.frame = m.frame + 1;
			fprintf('\rSaving frame: %d',m.frame); fflush(stdout);
		end
	else
		drawnow();
	end
end