function endMovieCapture(m,movieName = 'octaveMovie.avi')
	if inoctave()
		system(['ffmpeg -sameq -i plot-output/%05d.png -y ',m.loc,'/',movieName]);
		clf;
		fprintf('\rFinished saving movie\n');
	end
end