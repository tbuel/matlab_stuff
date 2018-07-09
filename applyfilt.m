function [varargout] =applyfilt(data,fs,fcut,p)
% function [varargout] =applyfilt(data,fs,fcut,p)
%
% data = input data
% fs = sample rate (if not specified, assumed 1sample/s)
% fcut = cutoff frequency (if not specified, fcut = fs/10)
% p = number of poles (if not specified, p = 2)
% 
% X = applyfilt(data,fs,fcut,p)
% X is the filtered output based on fs
%
% applyfilt(data,fs,fcut,p)
% by itself will plot all channels...
%
% Modular Script to Apply Filter
% SA Tbuel 9/30/2015

	% Input Checks
	nargoutchk(0,1);
	if p == [],
		p = 2;
	end
	if fcut == [],
		fcut = fs/10;
	end
	
	% Apply Filter per Column
	Nyq = fs/2;
	[M,N] = size(data);
	% in case we're dealing with a one row vector
	if(M == 1 && N > 1),
		data = data';
		[M,N] = size(data);
	end
	
	for n = 1:N,
		X = data(:,n);
		[B,A]=butter(p,fcut/Nyq);
		Xout(:,n) = filtfilt(B,A,X);
	end
	
	% If no arguments, plot the data
	if nargout==0,
		t = (0:M-1)/fs;
		for n = 1:N,
			titles = ['Filtered Ch' num2str(n); 'Fcutoff= ' num2str(fcut)];
			figure(n)
			plot(t,Xout(:,n));
			xlabel('Time, s')
			ylabel('Magnitude')
			title(titles)
		end
	end
	
	% Output Arguments
	varargout{1} = Xout;
end