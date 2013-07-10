function options = readOptions(varargin)
	options.dummyVar = 0; %create the structure

    pos = 1;
    while pos <= nargin
        field = varargin{pos};
        if ischar(field)
            options = setfield(options,field,varargin{pos+1});
            pos = pos + 2;
        elseif isstruct(field)
            options = combineStructs(options,field);
            pos = pos + 1;
        elseif iscell(field)
            options = combineStructs(options,readOptions(field{:}));
            pos = pos + 1;
        end
    end
end

function a = combineStructs(a,b)
    names = fieldnames(b);
    for i = 1:length(names)
        a = setfield(a,names{i}, getfield(b,names{i}));
    end
end