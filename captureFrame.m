function m = captureFrame(m)
	if inoctave()
		print(sprintf('%s/%05d.png',m.loc,m.frame));
		m.frame = m.frame + 1;
		fprintf('\rSaving frame: %d',m.frame); fflush(stdout);
	end
end