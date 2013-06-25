%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct the t values 
%
% build a mesh for sampling at constant t intervals
% there are more x-samples near the shoreline, this captures that

function [x_lin, eta_lin, u_lin, t_lin] = toConstantTime(x2, t2, eta2, u2, timeSamples)
    x_sample = mean(x2,2);
    x_sample = x_sample - x_sample(end);
    x_sample = x_sample ./ x_sample(1);

    % each t sample is the average of t(lambda = i)
    if exist('timeSamples')
        t_lin = timeSamples;
    else
        t_lin = mean(t2);
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

    % sample the wave height (eta) at the correct t-values
    eta_lin = griddata(v(x2), v(t2), v(eta2), x_lin, t_lin); % 'cubic'
    if exist('u2')
        u_lin = griddata(v(x2), v(t2), v(u2), x_lin, t_lin); 
    end
end