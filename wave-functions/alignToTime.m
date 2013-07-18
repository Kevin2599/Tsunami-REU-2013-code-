%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x_mesh, t_samples, varargout] = toConstantTime(x2, t2, timeSamples, [options,] varargin)
% function sout = toConstantTime(dataStruct, options)
%
% Correct the t values 
%
% timeSamples = {'mean','linear', [your values]}
% there are more x-samples near the shoreline, this captures that

function [x_mesh, t_samples, varargout] = alignToTime(x2, t2, timeSamples, varargin)

    options.dummyVar = [];
    if isstruct(varargin{1})
        options = varargin{1};
        varargin = varargin(2:end);
    end

  % specify the time samples
    if ~ischar(timeSamples)
        t_samples = timeSamples;
    elseif strcmp(timeSamples,'mean') % each t sample is the average of t(lambda = i)
        t_samples = mean(t2);
    elseif strcmp(timeSamples,'linear')
        t_samples = linspace( max(min(t2')), min(max(t2')), size(t2,2) );
    else
        error(['option ''' timeSamples ''' is not supported']);
    end

    x_mesh = interpOver(t2,x2,t_samples);
    for i = 1:(nargout-2)
        varargout{i} = interpOver(t2,varargin{i},t_samples);
    end

    xSamples = readOption(options,'toConstantTime_xSamples','default');
    if ~ischar(xSamples)
        % FIXME
    end
end

function YI = interpOver(X, Y, XI)
    YI = zeros(size(X,1),length(XI));
    for m=1:size(X,1)
        YI(m,:) = interp1(X(m,:),Y(m,:),XI);
    end
end