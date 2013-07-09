% function plotWave(x,height,bathymetry)

function plotWave(x,z,t,bath,varargin)
	options = readOptions(varargin);
    getOption = @(name,defaultValue) readOption(options,name,defaultValue);

    if getOption('plotTime',false)
        timePlot = makePlot('saveMovie',getOption('saveTimePlot',false),'movieName','waterOutline.avi', 'movieLocation','octaveMovies', ...
                                'subplot',[4 1],'hold','off');
        % time3DPlot = makePlot('hold', 'on');

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
        x_axis   = [min(min(x))  , max(max(x))+10];
        eta_axis = [min(min(z)), max(max(z)) ];

        [outlineInitial_x outlineInitial_y] = topViewOfWater(bath,[0; min(min(x))],[0; 0]);

    %% top view
        % switchToPlot(time3DPlot);
        % [bath_x bath_y] = meshgrid(x_axis,[bath.left; bath.right(end:-1:1)]);
        % surf(bath_x,bath_y,[bath.height bath.height ; bath.height(end:-1:1) bath.height(end:-1:1)] + ones(4,1) * (x_axis*bath.slope),'EdgeColor','none','LineStyle','none');

        for i=1:150 %length(t)

        %% top view
            switchToPlot(timePlot,1:2);

            [outline_x outline_y outline_dy] = topViewOfWater(bath,x(:,i),z(:,i), getOption('waveOutlineMagnification',1));

            plot(outlineInitial_x,outlineInitial_y,'k');
            plot(outline_x,outline_y,'r');

            xlim(x_axis);
            ylabel('Y');
            title(['top view, t=', num2str(t(i))]);

        %% displacement
            switchToPlot(timePlot,3);

            plot(outline_x(1:end/2),outline_dy);
            xlim(x_axis);
%            ylim(2.2 * eta_axis / DJN_slopes);
            ylabel('edge runup');
            grid on

        %% side view
            switchToPlot(timePlot,4);
            % plot(x2(:,i),eta2(:,i)', 'b')
            plot(x(:,i),z(:,i), 'r');

            plot(x_axis , bath.slope*x_axis);
            plot(0,0,'^b');

            axis([x_axis eta_axis]);
            leg=legend('lambda','real t');
            set(leg,'Location','northwest')
            xlabel(['X']); ylabel(['Z']);
            

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
            timePlot = drawPlot(timePlot);

        end
        finishPlot(timePlot);
    end
end