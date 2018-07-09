function [varargout] = routh(a)
% function [varargout] = routh(a)
%   Routh Stability Criterion
    if ~isvector(a)
        error('a must be a vector');
    end
    if iscolumn(a),
        a = a';
    end
    N = length(a);
    if a(1) < 0,
        a = -a;
    end
    M = zeros(N,N);
    cols = ceil(N/2);
    rts = roots(a);
    % First 2 Rows are Coefficients 
    M(1,1:cols) = a(1:2:N);
    if mod(N,2)
        M(2,1:cols-1) = a(2:2:N);
    else
        M(2,1:cols) = a(2:2:N);
    end
    nnp = 0; %nonnegative poles
    % Find rest of rows 
    for i = 3:N,
        d = M(i-1,1);
        for j = 1:cols,
            mdet = (M(i-1,1)*M(i-2,j+1) - M(i-2,1)*M(i-1,j+1));
            x = mdet / d;
            M(i,j) = x;
            if j == 1,
                if sign(M(i,j)) ~= sign(M(i-1,j))
                    nnp = nnp + 1; % Sign changes => more RHP poles
                end
            end
        end
    end
    
    if nargout == 0,
        fprintf('\nRouth-Hurwitz Array:\n');
        disp(M); % This actually looks way better
%         for i = 1:N,
%             fprintf('%d ',M(i,:));
%             fprintf('\n');
%         end
        fprintf('\nRHP Poles: %d\nSystem is ',nnp);
        if (nnp == 0),
            fprintf('stable\n'); 
        else
            fprintf('unstable\n'); 
        end
        fprintf('\nRoots:\n');
        disp(rts);
    else
        varargout{1} = M;
        varargout{2} = nnp;
        varargout{3} = rts;
    end

end

