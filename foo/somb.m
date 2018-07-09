function Z = somb(x,y)
% function z = somb(x,y)
%   somb(x) function.
%   x and y are 1xN and 1xM arrays

    
    if nargin == 0,
        N = 127;    
        x = linspace(-10,10,N);
        y = x;
    elseif nargin == 1,
        y = x;
    end
    [X,Y] = meshgrid(x,y);
    R = sqrt(X.^2+Y.^2);
    Z = besselj(1,R)./R;
    if nargout == 0
        figure
        surf(X,Y,Z)
    end
end

