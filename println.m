function println(text)
	if inoctave()
		fprintf([text '\n']);
		fflush(stdout);
	else
		disp(text);
	end
end