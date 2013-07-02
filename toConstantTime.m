%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x_lin, t_lin, varargout] = toConstantTime(x2, t2, timeSampleType, varargin)
%
% Correct the t values 
%
% timeSampleType = {'mean','linear'}
% build a mesh for sampling at constant t intervals
% there are more x-samples near the shoreline, this captures that

function [x_lin, t_lin, varargout] = toConstantTime(x2, t2, timeSampleType, varargin)
    x_sample = mean(x2,2);
    x_sample = x_sample - x_sample(end);
    x_sample = x_sample ./ x_sample(1);

    % each t sample is the average of t(lambda = i)
    if strcmp(timeSampleType,'mean')
        t_lin = mean(t2);
    elseif strcmp(timeSampleType,'linear')
        t_lin = linspace( max(min(t2')), min(max(t2')), length(t2(1,:)) );
    else
        error(['option ''' timeSampleType ''' is not supported']);
    end

    % moving shoreline for x
    x_max = interp1(t2(  1,:), x2(  1,:), t_lin);
    x_min = interp1(t2(end,:), x2(end,:), t_lin);

    % create the mesh
    t_lin = ones(size(x_sample)) * t_lin;
    x_lin = x_sample * x_max + (1-x_sample) * x_min;
    mesh(t_lin, x_lin, ones(size(x_lin)));

    % vectorize a matrix
    v = @(mat) reshape(mat(),1,[]);
    x2 = v(x2); t2 = v(t2);

    % sample the vars at the correct t-values
    for i = 1:(nargout-2)
        varargout{i} = griddata(x2, t2, v(varargin{i}), x_lin, t_lin); % 'cubic'
    end
end