% function [varargout] =ride_quality([X Y Z])
% Ride Quality (1/3 Octave RMS) Parameters
%
% SA Tbuel 5/23/17
function [out] = ride_quality(data,fs)

	[M,N] = size(data);
	g = 9.81; % 1g = 9.81m/s^2
	if N ~= 3,
		error('\nData must be a M x 3 matrix with Longitudinal, Lateral, and Vertical accelerations respectively\n');
	else
		load('f_low_3.dat'); % Lower Cutoff Frequencies
		load('f_high_3.dat'); % Higher Lower Cutoff Frequencies 
		Nyq = fs/2; % For next line only
		f_3lim = find(f_high_3 <= Nyq,1,'last');
		f_low = f_low_3((1:f_3lim));
		f_high = f_high_3((1:f_3lim));
		BF = length(f_low);

		%% 1/3 Octave Ride Quality
		for n = 1:3,
			for bf = 1:BF,
				Wn = [f_low(bf),f_high(bf)]/Nyq;
				[B,A]=butter(2,Wn);
				octn = filtfilt(B,A,g*data(:,n));
				octrms(bf,n) = rms(octn);	
			end
		end
		out = octrms;
	end
				
	
