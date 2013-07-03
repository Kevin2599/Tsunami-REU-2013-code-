function finishPlot(m)
	if m.saveMovie
		if inoctave()
			system(['ffmpeg -sameq -i ',m.frameLoc,'/%05d.png -y ',m.loc,'/',m.movieName]);
			system(['rm -rf ' m.frameLoc]);
			fprintf('\rFinished saving movie\n');
			graphics_toolkit fltk
		end
	end
	if readOption(m,'closeOnFinish',false)
		close(m.plotHandle);
	end
end