function R = rootsMat(C)
	for i=1:size(C,2)
		R(:,i) = roots(C(:,i));
	end
end
