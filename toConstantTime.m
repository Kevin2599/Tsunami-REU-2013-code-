%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x_mesh, t_samples, varargout] = toConstantTime(x2, t2, timeSamples, varargin)
%
% Correct the t values 
%
% timeSamples = {'mean','linear', [your values]}
% build a mesh for sampling at constant t intervals
% there are more x-samples near the shoreline, this captures that

function [x_mesh, t_samples, varargout] = toConstantTime(x2, t2, timeSamples, varargin)
    x_sample = mean(x2,2);
    x_sample = x_sample - x_sample(end);
    x_sample = x_sample ./ x_sample(1);

  % each t sample is the average of t(lambda = i)
    if ~ischar(timeSamples)
        t_samples = timeSamples;
    elseif strcmp(timeSamples,'mean')
        t_samples = mean(t2);
    elseif strcmp(timeSamples,'linear')
        t_samples = linspace( max(min(t2')), min(max(t2')), length(t2(1,:)) );
    else
        error(['option ''' timeSamples ''' is not supported']);
    end

  % moving shoreline for x
    x_max = interp1(t2(  1,:), x2(  1,:), t_samples);
    x_min = interp1(t2(end,:), x2(end,:), t_samples);

  % create the mesh
    t_mesh = ones(size(x_sample)) * t_samples;
    x_mesh = x_sample * x_max + (1-x_sample) * x_min;

  % plot the mesh
    % mesh(t_mesh, x_mesh, ones(size(x_mesh)));

  % vectorize a matrix
    v = @(mat) reshape(mat(),1,[]);
    x2 = v(x2); t2 = v(t2);

  % sample the vars at the correct t-values
    for i = 1:(nargout-2)
        varargout{i} = griddata(x2, t2, v(varargin{i}), x_mesh, t_mesh); % 'cubic'
    end
end