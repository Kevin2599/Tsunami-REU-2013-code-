function sout = applyFunToStruct(f,s)
	names = fieldnames(s);
	sout = [];
    for i = 1:length(names)
        sout = setfield(sout,names{i}, f(getfield(s,names{i})));
    end
end