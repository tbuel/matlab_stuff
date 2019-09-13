function [varargout] = pstrain(varargin,opts)
% function [varargout] = pstrain(varargin,opts)
% Computes principal strains/stresses from rosette data
% Currently computes p-strains rectangular rosette only
% opts = [] : computes strain
% opts = 'stress' : computes stress

    rad2 = sqrt(2);
    E = 29e9; % elastic modulus
	nu = 0.26; % poisson's ratio
	
    if nargin == 1,
        [M,N] = size(varargin);
        if N ~= 3,
            error('Invalid input matrix size. \nEither input 3 distinct vectors, or a single Mx3 matrix');
        end
        data = varargin{1};
        strain0 = data(:,1);
        strain45 = data(:,2);
        strain90 = data(:,3);
    elseif nargin == 3 || nargin == 4,
		% if nargin == 3,
			strain0 = varargin{1};
			strain45 = varargin{2};
			strain90 = varargin{3};
			M = length(strain0);
		% end
    else
        M = 2^10;
        fs = 1e3;
        t = (0:M-1)/fs;
        strain0 = 50*cos(2*pi*20.*t) + 5*rand([1,M]);
        strain45 = 10*cos(2*pi*10.*t + pi/4) + 5*rand([1,M]);
        strain90 = 50*cos(2*pi*20.*t + pi/2) + 5*rand([1,M]);
    end
    
	s_determinant = sqrt((strain0 - strain45).^2 + (strain45 - strain90).^2);
	phi = 0.5*atan((strain0-2*strain45+strain90)./(strain0-strain90)); % Angle
    if opts == [],
		% Strain output - Default
		strain1 = (strain0 + strain90)/2 + s_determinant/rad2;
		strain2 = (strain0 + strain90)/2 - s_determinant/rad2;
	elseif strcmp(opts,'stress'),
		% Stress output - variable names same for easy
		strain1 = 0.5*E*((strain0 + strain90)/(1-nu) + rad2*s_determinant/(1+nu));
		strain2 = 0.5*E*((strain0 + strain90)/(1-nu) - rad2*s_determinant/(1+nu));
	else 
		error('Invalid option');
	end
    
    % If no output, plot Principal Strains
    if nargout == 0,
%         figure('units','normalized','outerposition',[0 0 1 1]),
        xlab = 'Samples';
        if nargin == 0,
            t = 1:M;
            xlab = 'Time, s';    
        end
        figure()
        subplot(2,2,1)
        stitle = 'Input Strains - Time History';
        plot(t,[strain0; strain45;strain90]) %,'YTick',ticks)
        xlabel(xlab) %,'Color','black','FontSize',label_size,'FontWeight','bold'
        ylabel('Strain, uS')
        legend(['Strain 0' char(176)],['Strain 45' char(176)],['Strain 90' char(176)])
        title(stitle,'FontSize',16)
        set(gca,'FontSize',12);
        grid on
        ysym(gca);
        
        
        subplot(2,2,2)
        stitle = 'Principal Strains - Time History';
        plot(t,[strain1; strain2]) %,'YTick',ticks)
        xlabel(xlab) %,'Color','black','FontSize',label_size,'FontWeight','bold'
        ylabel('Strain, uS')
        legend('Principal Strain 1','Principal Strain 2')
        title(stitle,'FontSize',16)
        set(gca,'FontSize',12);
        grid on
        ysym(gca);
        
        subplot(2,2,3)
        stitle = 'Principal Strains - Angle';
        plot(t, phi) %,'YTick',ticks)
        xlabel(xlab) %,'Color','black','FontSize',label_size,'FontWeight','bold'
        ylabel('Angle, rad')
        title(stitle,'FontSize',16)
        set(gca,'FontSize',12);
        ysym(gca);
        grid on
        
        
        subplot(2,2,4)
        stitle = 'Principal Strains - XY Plot';
        plot(strain1, strain2) %,'YTick',ticks)
        xlabel('Principal Strain 1, uS') %,'Color','black','FontSize',label_size,'FontWeight','bold'
        ylabel('Principal Strain 2, uS')
        title(stitle,'FontSize',16)
        set(gca,'FontSize',12);
        grid on
	elseif nargout == 1,
		out(:,1) = strain1;
		out(:,2) = strain2;
		out(:,3) = phi;
		varargout{1} = out;
	elseif nargout == 3,
		varargout{1} = strain1;
		varargout{2} = strain2;
		varargout{3} = phi;
    end
end

