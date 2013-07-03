function switchToPlot(m)
	if gcf() ~= m.plotHandle
		figure(m.plotHandle)
	end
end