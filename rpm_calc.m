%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPM Calculator
% Code Taken from NS_MX_Torq_Accel.m
% Tanner 02/24/2015
% Requires Signal Processing Toolbox
%
% To Do:
% -Correlate variables appropriately
% -Multiple Column Compatibility
% -Function Call Error Handling
% -Optimize
% -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rpmout = rpm_calc(X,fs,offset,scale,maxrpm,min_pks)
	[M,N] = size(X);
	
	maxrpm = 3000; %max RPM
	min_pkd = fs*60/maxrpm; % This is for the RPM routine... 10 samples
    
	rpmtmp = data(:,14);
    rpmtmp = mean(rpmtmp)-rpmtmp; % 01/12/2015
    min_pks = 0.3; % 01/12/2015
	
    [~,locs1] = findpeaks(rpmtmp,'MINPEAKHEIGHT',min_pks);
    [~,locs2] = findpeaks(rpmtmp,'MINPEAKDISTANCE',min_pkd);
	
    e = 1;
    for c = 1:length(locs1)
        for d = 1:length(locs2)
            if (locs1(c) == locs2(d)),
                locs(e) = locs2(d);
                e = e+1;
                break;
            end
        end
    end
    clear locs1
    clear locs2
    
    rpmact = zeros(length(locs)-1,2);
    rpmact(:,1) = (locs((2:length(locs)))-1)/fs;
    rpmact((2:length(rpmact)),2) = 60./diff(rpmact(:,1));
    %filerpm = strcat(filename,'_rpm.txt')
    %save(filerpm,'-ascii','-tabs','rpmact');
    clear rpmact
    clear filerpm
    clear rpmtmp %1/5/2015
    clear min_pks
    clear mjr_pks
    clear locs %1/5/2015