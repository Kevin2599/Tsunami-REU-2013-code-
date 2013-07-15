% function plotWave(x,height,bathymetry)

function plotWave(data,bath,varargin)
	options = readOptions(varargin);
    getOption = @(name,defaultValue) readOption(options,name,defaultValue);

    x = data.x; t = data.t; z = data.eta;

    x_axis  = [min(min(x)), max(max(x)) ];
    x_axis(2) = x_axis(2) + .05 * (x_axis(2) - x_axis(1));
    z_axis  = [min(min(z)), max(max(z)) ];

    if getOption('plotBathymetry',false)
        plotWindow = makePlot('timePlot','saveMovie',getOption('saveTimePlot',false),'movieName','waterOutline.avi', 'movieLocation','octaveMovies', ...
                            'subplot',[4 1],'hold','off');
        bathPlot = makeSubplot(plotWindow,1:2);
        dyPlot   = makeSubplot(plotWindow,3);
        zPlot    = makeSubplot(plotWindow,4);

        % max_slope = max(diff(bath.right - bath.left) ./ diff(bath.height));
        % disp_axis = z_axis * max_slope;

        x_samples = x(:,1) ./ (-1*min(min(x)));
        [outlineInitial_x outlineInitial_y] = topViewOfWater(bath,x_samples,0*x_samples);

        % [bath_x bath_y] = meshgrid(x_axis,[bath.left; bath.right(end:-1:1)]);
        % bath_z = [bath.height bath.height ; bath.height(end:-1:1) bath.height(end:-1:1)] + ones(4,1) * (x_axis*bath.slope);
        % y_axis = [min(min(bath_y)) max(max(bath_y))];
    else
        plotWindow = makePlot('timePlot','saveMovie',getOption('saveTimePlot',false),'movieName','waterOutline.avi', 'movieLocation','octaveMovies', ...
                            'hold','off');
        zPlot = plotWindow;
    end

%% line plot of x,t,eta
    if ~isvector(t) && getOption('plotTime_viewSamples',false)
        time3DPlot = makePlot('time3DPlot','hold', 'on');
        switchToPlot(time3DPlot);
        for i=1:length(x(1,:))
            plot3(t(:,i)', x(:,i)', z(:,i)');
            time3DPlot = drawPlot(time3DPlot);
        end
        surf(x,t,z,'EdgeColor','none','LineStyle','none');
    end

    timeValues = 1:size(x,2);
    if getOption('limitPlotT',false)
        timeValues = options.limitPlotT;
    end

    for i=timeValues
    %% top view
        if exist('bathPlot')
            bathPlot = switchToPlot(bathPlot);

            [outline_x outline_y outline_dy] = topViewOfWater(bath,x(:,i),z(:,i), getOption('waveOutlineMagnification',1));

            plot(outlineInitial_x,outlineInitial_y,'k');
            plot(outline_x,outline_y,'r');

            xlim(x_axis);
            ylabel('Y');
            title(['top view, t=', num2str(t(i))]);
        end

    %% displacement
        if exist('dyPlot')
            dyPlot = switchToPlot(dyPlot);

            plot(outline_x(1:end/2),outline_dy);

            xlim(x_axis);
            % ylim(disp_axis);
            ylabel('edge runup');
            grid on
        end

    %% side view
        if exist('zPlot')
            zPlot = switchToPlot(zPlot);
            % plot(x2(:,i),eta2(:,i)', 'b')
            plot(x(:,i),z(:,i));%, '.r');

            plot(x_axis , bath.slope*x_axis,'r');
            plot(0,0,'^b');

            axis([x_axis z_axis]);
            % set(leg,'Location','northwest')
            xlabel(['X']); ylabel(['Z']);
            if ~isvector(t)
                title(['t = ' num2str(t(2,i))]);
            end
        end

        plotWindow = drawPlot(plotWindow);

    %% 3D plot
        if getOption('timePlot3D',false)
            time3DPlot = switchToPlot(time3DPlot);

            surf(bath_x,bath_y, bath_z,   2*ones(size(bath_x)),'EdgeColor','none','LineStyle','none');
            surf([1 ; 1] * x(:,i)', y_axis' * ones(size(x(:,i)')), [1 ; 1 ] * z(:,i)', ...
                    ones(2,length(x(:,i))),'EdgeColor','none','LineStyle','none')

            time3DPlot = drawPlot(time3DPlot);
        end
        
    end
    finishPlot(plotWindow);
end


%% OLD PLOTTING FUNCTIONS


        % figure(4); hold off
        % semilogy(x(:,i),J(:,i));
        % axis([x_axis  min(min(J)) max(max(J))]);
        % title('Jacobian');

        % xlabel(['Jacobian=',num2str(min(J(:,i)))]);
        % if min(J(:,i))<0
        %     pause(.1)
        % end
        % pause(framerate);
        % figure(3,'visible','on')

% %     % Plot to look for global error and information
    % slope=zeros(1,length(lambda));
    % breakc=slope;
    % for j=1:length(lambda)
    %     slope = diff(eta2(:,j)) ./ diff(x2(:,j));
    %     % for i=1:length(eta2(:,1))-1
    %     %     slope(i)=(eta2(i+1,j)-eta2(i,j))/(x2(i+1,j)-x2(i,j));
    %     % end
    %     breakc(j)=max(slope(:));
    % end
    % for i=1:length(lambda)
    %     % %         if ((breakc(i)>=1/2*bath.slope)||(i==brokeat))
    %     % %             println('BROKE...')
    %     % %             if found
    %     % %                 println('Numerical')
    %     % %             end
    %     % %             %break
    %     % %         end
    %     index1=(J(:,i)>=0);
    %     index2=~index1;
    %     plot(sigma(index1), eta2(index1,i), '.r')
    %     hold on
    %     plot(sigma(index2), eta2(index2,i), '.b')
    %     plot(sigma, J(:,i), '-k')
    %     hold off
    %     axis([0 300 min(min(eta2)) max(max(eta2))])
    %     leg=legend('Aprox solution');
    %     % set(leg,'Location','Best');
    %     xlabel('Sigma')
    %     title(num2str(t2(2,i)))
    %     pause(0.01)
    % end
