function [varargout] = mup_speed_convert(data,fs,DiaW,Nteeth,Threshold)
% function [varargout] = mup_speed_convert(data,fs,DiaW,Nteeth,Threshold)
%
% data = input data, in mV
% fs = sample rate, Hz
% DiaW = Wheel Diameter, inches
% Nteeth = number of teeth on gear
% Threshold = threshold for checking input signal, mV (default: 1% of Max)
% 
% Speed = mup_speed_convert(data,fs,DiaW,Nteeth,Threshold)
% 'Speed' is the wheel's speed in mph
%
% mup_speed_convert(data,fs,DiaW,Nteeth,Threshold)
% by itself will plot all channels...
%
% Modular Script to convert MUP speed signal to Speed, mph
% SA Tbuel 6/18/2018

	% Input Checks
	nargoutchk(0,2);
    
	if Threshold == [],
		Threshold = 0.01*max(data);
	end

    M = length(data);
    freq2mph = 0.1785 * DiaW / Nteeth; % Conversion Factor
    
    tspeed = 0:floor(M/fs); % time vector
    speednew = zeros(length(tspeed),1); % new speed in mph
    
    i = 1; % increment
    zerox = 0; % number of zero crossings

    for m = 1:fs:M,
        if m + fs - 1 > M,
            last = M;
        else
            last = m+fs-1;
        end
        if rms(data(m:last)) > Threshold,
            if m == 1,
                continue;
            end
            for mm = m:last,
                if sign(data(mm)) ~= sign(data(mm-1)),
                    zerox = zerox + 1;
                end
            end
        end
        freq1s = zerox/2;
        speednew(i) = freq1s*freq2mph;
        i = i + 1;
        zerox = 0;
    end        
    	
	% If no arguments, plot the data
	if nargout==0,
        figure(n)
        plot(tspeed,speednew);
        xlabel('Time, s')
        ylabel('Speed, mph')
        title('Speed of MUP Sensed Wheel')
    % Output Variables    
    elseif nargout == 1,
        varargout{1} = speednew;
    elseif nargout == 2,
        varargout{1} = speednew;
        varargout{2} = tspeed';
    end
        
end

