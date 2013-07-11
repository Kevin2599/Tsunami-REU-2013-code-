%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x_mesh, t_samples, varargout] = toConstantTime(x2, t2, timeSamples, varargin)
%
% Correct the t values 
%
% timeSamples = {'mean','linear', [your values]}
% build a mesh for sampling at constant t intervals
% there are more x-samples near the shoreline, this captures that

function [x_mesh, t_samples, varargout] = toConstantTime(x2, t2, timeSamples, varargin)
    options.dummyVar = [];
    if isstruct(varargin{1})
        options = varargin{1};
        varargin = varargin(2:end);
    end

    xSamples = readOption(options,'toConstantTime_xSamples','default');
    if strcmp(xSamples,'default')
        x_sample = mean(x2,2);
    elseif ~ischar(xSamples)
        x_sample = xSamples
    end

    x_sample = x_sample - x_sample(end);
    x_sample = x_sample ./ x_sample(1);

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
        varargout{i} = griddata(x2, t2, v(varargin{i}), x_mesh, t_mesh); % other option, 'cubic'

        if ~inoctave()
            [e x t] = DJN_bluh(x2,t2,varargin{i});
            figure(1);
            mesh(x,t,e);
            figure(2);
            mesh(x_mesh,t_mesh,varargout{i});
        end
    end
end

function [e x t] = DJN_bluh(x2,t2,eta2)
    x21=x2(2:end,:);
    t21=t2(2:end,:);
    eta21=eta2(2:end,:);
    F=TriScatteredInterp(x21(:), t21(:), eta21(:));

    min_x=-500;         max_t=500;
    dx=0.1;
    dt=0.1;

    [x,t]=meshgrid(min_x:dx:max(x21(:))*1.1, 0:dt:max_t);
    e=F(x,t);

    index=(x>-50);
    for i=find(index==1)
        z=interp1(t21(1,:), x21(1,:), t(i));
        if(x(i)>z)
            e(i)=NaN;
        end
    end

    line_x=[x21(1,:) min_x min_x  x21(1,1)];
    line_y=[t21(1,:) max_t     0  t21(1,1)];

    index=inpolygon(x, t, line_x, line_y);
    e(~index)=NaN;
end