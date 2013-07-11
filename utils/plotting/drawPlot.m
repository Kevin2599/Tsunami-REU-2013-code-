function m = drawPlot(m)
	if m.saveMovie
		if inoctave()
			fileName = sprintf('%s/%05d.png',m.frameLoc,m.frame);
			print(fileName);
			m.frame = m.frame + 1;
			fprintf('\rSaving frame: %d',m.frame); fflush(stdout);
		else
			frame = getframe();
			writeVideo(m.writer,frame);
		end
	else
		refresh();
	end
end