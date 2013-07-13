function m=switchToPlot(m,subplotIndex)
	if m.plotHandle == -1
		global plotIDs
		if ~exist('plotIDs')
			plotIDs.dummy__ = 0;
		end

		if isfield(plotIDs,m.plotID)
			figure(getfield(plotIDs,m.plotID));%,'visible',readOption(m,'visible','on')); clf('reset');
		else
			figure();%'visible',readOption(m,'visible','on')); clf('reset');
			plotIDs = setfield(plotIDs,m.plotID,gcf());
		end
		m.plotHandle = gcf();
		if isfield(m,'subplot')
			m.subplotDim = m.subplot;
			subplot(m.subplotDim(1), m.subplotDim(2), 1)
		end
	elseif gcf() ~= m.plotHandle
		figure(m.plotHandle)
	end
	if exist('subplotIndex')
		subplot(m.subplotDim(1),m.subplotDim(2),subplotIndex);
	end
	if isfield(m,'hold')
		hold(m.hold)
		x = xlim();
		y = ylim();
		plot(mean(x),mean(y),'.','MarkerSize',.001);
		hold on;
	end
end