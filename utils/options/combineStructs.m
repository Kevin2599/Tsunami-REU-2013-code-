function a = combineStructs(a,b,varargin)
	if nargin > 2
		b = combineStructs(b,varargin{:});
	end
    names = fieldnames(b);
    for i = 1:length(names)
        a = setfield(a,names{i}, getfield(b,names{i}));
    end
end