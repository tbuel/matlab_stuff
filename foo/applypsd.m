%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [varargout] =applypsd(data,fs,NFFT)
% data = input data
% fs = sample rate [samples/s]
% NFFT = number of points of which to do FFT (default = 2^(nextpow2(nyquist))
%
% Modular Script to do PSD based on PWELCH estimate
% SA Tbuel 9/30/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] =applypsd(data,fs,NFFT)
	nargoutchk(0,2); % Argument Check
	[M,N] = size(data);
	if NFFT == [],
		NFFT = 2^nextpow2(10*fs); % If not defined, use 2^nextpow2(Nyquist)
        disp(NFFT)
	end
	% 1-row vector
	if M == 1 && N > 1,
		data = data';
		[M,N] = size(data);
	end
	%% Power Spectral Density
	for n = 1:N,
		[Pxx,f] = pwelch(data(:,n),[],[],NFFT,fs);
		if nargout == 2,
			PSDout(:,n) = Pxx;
		else
			if n == 1,
				PSDout(:,1) = f;
			end
			PSDout(:,n+1) = Pxx;
		end	
	end
		
	% Output Variables
    if nargout == 0,
        % plot if no outputs
        for n = 1:N
			figure(n)
			semilogx(f,10*log10(Pxx(:,n)));
			xlabel('Frequency, Hz')  
			ylabel('PSD, dB/Hz')  
			titles = ['Power Spectral Density of Ch' num2str(n)];
			title(titles)
            grid on
        end	
    elseif nargout == 1,
        varargout{1} = PSDout;
    elseif nargout == 2,
        varargout{1} = f;
        varargout{2} = PSDout;
    end
end
