% function plotWave(x,height,bathymetry)

function plotWave(x,z,t,bath,varargin)
	options = readOptions(varargin);
    getOption = @(name,defaultValue) readOption(options,name,defaultValue);

    if isvector(t) %% x(t), z(t)
        timePlot = makePlot('timePlot','saveMovie',getOption('saveTimePlot',false),'movieName','waterOutline.avi', 'movieLocation','octaveMovies', ...
                                'subplot',[4 1],'hold','off');
        time3DPlot = makePlot('time3DPlot','hold', 'on');

    %% line plot of x,t,eta
    %% warning, this requires a matrix for t
        if getOption('plotTime_viewSamples',false)
            figure(2); clf
            hold on
            for i=1:length(x2(1,:))
                plot3(t2(:,i)', x2(:,i)', eta2(:,i)');
                drawnow();
            end
            surf(x,t,z,'EdgeColor','none','LineStyle','none');
        end

    %% Bathymetry
        x_axis    = [min(min(x))  , max(max(x))+10];
        eta_axis  = [min(min(z)), max(max(z)) ];

        max_slope = max(diff(bath.right - bath.left) ./ diff(bath.height));
        disp_axis = eta_axis * max_slope;

        [outlineInitial_x outlineInitial_y] = topViewOfWater(bath,[0; min(min(x))],[0; 0]);

        [bath_x bath_y] = meshgrid(x_axis,[bath.left; bath.right(end:-1:1)]);
        bath_z = [bath.height bath.height ; bath.height(end:-1:1) bath.height(end:-1:1)] + ones(4,1) * (x_axis*bath.slope);
        y_axis = [min(min(bath_y)) max(max(bath_y))];

        for i=1:length(t)

        %% top view
            timePlot = switchToPlot(timePlot,1:2);

            [outline_x outline_y outline_dy] = topViewOfWater(bath,x(:,i),z(:,i), getOption('waveOutlineMagnification',1));

            plot(outlineInitial_x,outlineInitial_y,'k');
            plot(outline_x,outline_y,'r');

            xlim(x_axis);
            ylabel('Y');
            title(['top view, t=', num2str(t(i))]);

        %% displacement
            timePlot = switchToPlot(timePlot,3);

            plot(outline_x(1:end/2),outline_dy);

            xlim(x_axis);
            ylim(disp_axis);
            ylabel('edge runup');
            grid on

        %% side view
            timePlot = switchToPlot(timePlot,4);
            % plot(x2(:,i),eta2(:,i)', 'b')
            plot(x(:,i),z(:,i), 'r');

            plot(x_axis , bath.slope*x_axis);
            plot(0,0,'^b');

            axis([x_axis eta_axis]);
            leg=legend('real t');
            set(leg,'Location','northwest')
            xlabel(['X']); ylabel(['Z']);


            timePlot = drawPlot(timePlot);

        %% 3D plot
            if getOption('timePlot3D',false)
                time3DPlot = switchToPlot(time3DPlot);

                surf(bath_x,bath_y, bath_z,   2*ones(size(bath_x)),'EdgeColor','none','LineStyle','none');
                surf([1 ; 1] * x(:,i)', y_axis' * ones(size(x(:,i)')), [1 ; 1 ] * z(:,i)', ...
                        ones(2,length(x(:,i))),'EdgeColor','none','LineStyle','none')

                time3DPlot = drawPlot(time3DPlot);
            end
            

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
        end
        finishPlot(timePlot);


    else %%  x(s), z(s), t(s)
    	
        lambdaPlot = makePlot('lambdaPlot', 'hold','off');
        topPlot    = makePlot('lambdaPlotTop', 'hold','on');

        x_axis = [min(min(x)) max(max(x))];

        % Plot at the shore
        for i=1:length(x(1,:))
            if getOption('plotTopView',false)
            	topPlot = switchToPlot(topPlot);

                plot3(x(:,i)', t(:,i)', z(:,i)');

                topPlot = drawPlot(topPlot);
            end

            lambdaPlot = switchToPlot(lambdaPlot);

            plot(x(:,i),z(:,i), '.r');
            plot(x_axis, bath.slope*x_axis);
            plot(0,0,'^b');

            axis([-5, 1.5*max(max(x)), max(min(min(z)), -1), 1.5*max(max(z))]);

            xlabel('X'); ylabel('Z');
            title(['t = ' num2str(t(2,i))]);
            
            % results.snapshot{i}.x=x(:,i);
            % results.snapshot{i}.eta=z(:,i);
            % results.snapshot{i}.time=t2(2,i);
            % results.max_runup=max(max(eta2));
            % results.case=['case_',num2str(DJN_beachwidth),'m_',num2str(1/DJN_slopes),'_',num2str(bath.slope)];

            lambdaPlot = drawPlot(lambdaPlot);
        end
    end
end


%% OLD PLOTTING FUNCTIONS


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
