function switchToPlot(m,subplotIndex)
	if gcf() ~= m.plotHandle
		figure(m.plotHandle)
	end
	if isfield(m,'subplotDim')
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