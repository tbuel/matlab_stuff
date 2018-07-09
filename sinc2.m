function Z = sinc2(x,y)
% function z = sinc2(x,y)
%   Sinc^2(x) function.
%   x and y are 1xN and 1xM arrays

    if nargin == 0,
        N = 127;    
        x = linspace(-3,3,N);
        y = x;
    elseif nargin == 1,
        y = x;
    end
    [X,Y] = meshgrid(x,y);
    Z = sinc(sqrt(X.^2+Y.^2));
    if nargout == 0
        figure
        surf(X,Y,Z)
    end
end

