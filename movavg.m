% function out = movavg(data,window)
%
% Simple moving average - taken from Matlab example code for Filter
% data = input data
% window = number of points to average over
% Uses Filter function in matlab, with a = 1, and b(k) = 1/Window; k = [1,Window]
function out = movavg(data,window)
	% SA Tbuel 8/25/16
	b = ones(window,1)/window;
	out = filter(b,1,data);