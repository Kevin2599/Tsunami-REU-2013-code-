function options = readOptions(args)
	options.dummyVar = 0; %create the structure

    if length(args) ~= 0
        startPos = 1;
        if isstruct(args{1})
            options = args{1};
            startPos = 2;
        end
        for i=startPos:2:length(args)
            options = setfield(options,args{i},args{i+1});
        end
    end
end