% Peak 1-s RMS
% Tanner Buel 12/2014
%
% Find the Peak 1-s RMS of signal X with sampling rate fs
% function [vargout] = pk1srms(X,fs)
%
% 1. pkrms = pk1srms(X,fs)
%    'X' is a column-based data set with sampling rate 'fs',
%    then 'pkrms' will be the peak 1-s RMS 
%    = max(RMS) over fs number of points
% 
% 2. [inx, pkrms] = pk1srms(X,fs)
%    'inx' will be the first index (in time, s) where the
%    peak 1-s RMS was calculated 
%
% To Come... vargin with timestep option (default 10 steps)
function [varargout] = pk1srms(X,fs)
	nargoutchk(1,2); % Check Arguments
	[M, N] = size(X);
	inx = zeros(1,N);
	pkrms = zeros(1,N);
	% Peak 1-s RMS
	for n = 1:N,
		Xi = X(:,n);
		pktmp = 0;
		rmsinx = 0;
		for k = 1:10:M-fs,
			rmstmp = rms(Xi([k:k+fs]));
			if(pktmp < rmstmp)
				pktmp = rmstmp;
				rmsinx = k;
			end
		end
		pkrms(n) = pktmp;
		inx(n) = (rmsinx-1)/fs; % in Time - Time vector starts at 0
		clear Xi
		clear pktmp
		clear rmstmp
		clear rmsinx
	end
	% Output Arguments
	varargout{1} = pkrms;
	varargout{2} = [inx, pkrms];
