function [varargout] = gaussian(x,mu,sigma)
%function out = gaussian(x,mu,sigma)
%   Single gaussian pulse:
%   x = data vector
%   mu = mean
%   sigma = stddev
%   out = exp(-(x-mu).^2./(2*sigma^2) / (sqrt(pi.*2*sigma^2))
%
%Note: If mu or sigma are not specified, mu = 0, sigma = 1/sqrt(2*pi)

%% Checks
    if nargin < 3,
        sigma = 1/sqrt(2*pi);
        if nargin < 2,
            mu = 0;
        end
        if nargin == 0,
            x = linspace(-5*sigma,5*sigma,1024);
        end
    end
    if ~isvector(x),
        error('x must be an array.');
    end
    if ~isvector(sigma),
        error('sigma must be an array');
    end
    if iscolumn(x),
        x = x';
    end
    if iscolumn(sigma),
        sigma = sigma';
    end
    
%% Computation
    xvar2 = 2*sigma.^2;
    out = exp(-(x'-mu).^2./(xvar2))./(sqrt(pi.*xvar2));
    if nargout == 0,
        sig1 = sigma(1); % only show first sigma if applicable
        xtx = 5*sig1; % show +/-5 standard deviations
        figure
        plot(x,out(:,1));
        title({'Gaussian Pulse'; ...
              ['{\mu}=' num2str(mu) '; {\sigma}=' num2str(sig1) '; FWHM=' num2str(2.3548*sig1)] })
%         title({'Gaussian Pulse'; ...
%               ['$$\mu=' num2str(mu) '; \sigma=1/\sqrt{2\pi}; FWHM=' num2str(2.3548*sigma) '$$'] },'Interpreter','latex')
        xticks(-xtx:sig1:xtx)
        xticklabels({'-5{\sigma}','-4{\sigma}','-3{\sigma}','-2{\sigma}','-1{\sigma}','0'...
                     '1{\sigma}','2{\sigma}','3{\sigma}','4{\sigma}','5{\sigma}'})
        grid on
        xlim([-xtx xtx])
    elseif nargout == 1,
        varargout{1} = out;
    end
end

