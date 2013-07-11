function endMovieCapture(m,movieName)
    if ~exist('movieName')
        movieName = 'octaveMovie.avi';
    end
	if inoctave()
		system(['ffmpeg -sameq -i plot-output/%05d.png -y ',m.loc,'/',movieName]);
		clf;
		fprintf('\rFinished saving movie\n');
	end
end