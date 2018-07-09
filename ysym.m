function ysym(gca)
	% function ysym(gca)
    %    Centers plot on x-axis (y = 0)
	ylimits = get(gca,'Ylim');
	ymax = max(abs(ylimits));
	set(gca,'Ylim',[-ymax ymax]);
end