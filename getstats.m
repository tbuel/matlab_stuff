%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [varargout] =getstats(data)
%
% data = input data
% fs = sample rate (if not specified, assumed 1sample/s)
% fcut = cutoff frequency (if not specified, fcut = fs/10)
% p = number of poles (if not specified, p = 2)
% 
% X = getstats(data,fs)
% X(:,1) will be time vector starting at 0s
%
% [X t] = getstats(data,fs)
% should be used to not accept time
%
% Modular Script to Apply Filter
% SA Tbuel 9/30/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yout =getstats(data)
	nargoutchk(0,4); 	% Input Checks
	
	%% Max, Min, RMS
	% Only useful for first 13 channels but run it on all anyways
	ymin = min(data); %min
	ymax = max(data); %max
	yrms = rms(data); %RMS
	yavg = mean(data); % Arithmetic Mean
	% ypk1srms = pk1srms(data);
	
	yout = [ymin;ymax;yrms;yavg];
	% dataout = strcat(filename,'_stats.txt');
	% save  (dataout,'-ascii', '-tabs','yout');
	
	% Output Arguments
	varargout{1} = yout;
