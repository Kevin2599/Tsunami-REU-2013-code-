function value = readOption(options,name,defaultValue)
    if isfield(options,name)
        value = getfield(options,name);
    else
        value = defaultValue;
    end
end