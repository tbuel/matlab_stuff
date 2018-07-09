function [varargout] = kruskal(A,opt)
%function [B,wt] = kruskal(A,opt)
%kruskal() implements Kruskal's Algorithm for finidng 
%the minimum spanning tree. 
% 
%   A - Input Matrix
%   opt - Option for input matrix type: 
%         'Matrix' :(Default) A is assumed to be an NxN matrix with each
%         possible node combination accounted for. Here, N is the 
%         number of vertices/nodes in the network. The value at A(i,j) = WT. 
% 
%         'Vertices' : A is assumed to be an Nx3 matrix, with N as the 
%         number of connected nodes, or edges. The first column is the 
%         weight, the second column is the first node, i, and the third 
%         column is the second node, j. 
%   B - Minimum Spanning Tree Result
%   wt - New MST weight
% 
%If no output is given, this will output the MST tabular data and plot 
%the differences between networks.
    
% Initialization
    if nargin == 0,
        A = input('\nEnter Network:');
    end
    if nargin == 1,
        opt = 'Matrix';
    end
    
% Create list of Edges
    if strcmp(opt,'Matrix')
        [N,~]=size(A);
        E = numel(A(A>0))/2; % number of possible edges
        S = zeros(E,3);
        e = 1;
        for i = 1:N,
            for j = i+1:N,
                x = A(i,j);
                if x ~= 0,
                    S(e,:) = [x,i,j];
                    e = e + 1;
                end
            end
        end  
        wt0 = sum(sum(A))/2;
    elseif strcmp(opt,'Vertices')
        wt0 = sum(A(:,1))
        S = A
        N = max(max(A(:,2:3)))
    end
          
    S = sortrows(S,1); % Sort edges in "non-decreasing order of weight"
    
    % Algorithm  
    B = zeros(N-1,3);
    parent = 1:N;
    rank = zeros(1,N);
    i = 1;
    e = 1;
    while e <= N-1,
        edge = S(i,:); % Get smallest edge
        i = i + 1;
        x = sfind(parent,edge(2));
        y = sfind(parent,edge(3));
        if x ~= y,
            B(e,:) = edge;
            e = e + 1;
            [parent,rank] = subunion(parent,rank,x,y);
        end
    end
    wt = sum(B(:,1));
    
    % Output
    if nargout == 0,
        fprintf('\nMinimum Spanning Tree:\n\tWt\tNode1\tNode2\n');
        disp(B);
        fprintf('Input Tree Weight: %d\nMST Weight: %d\n',wt0,wt);
        
        figure
        phi = (0:N-1)*2*pi/N;
        rho = ones(1,N);
        % Spanning Tree
        for i = 1:length(S),
            phi1 = phi(S(i,2));
            phi2 = phi(S(i,3));
            polarplot([phi1 phi2],[1 1],':b')
            hold on
            avg = (phi1 + phi2)/2;
            if abs(phi2 - phi1) > pi
                avg = avg + pi;
                ravg = cos(avg-phi2);
            else
                ravg = cos(avg-phi1);
            end
            text(avg,1.1*ravg,num2str(S(i,1)))
            if i < N,
                polarplot([phi(B(i,2)) phi(B(i,3))],[1 1],'r')
            end
        end
        
        
        
%         for i = 1:N,
%             phi1 = phi(i);
%             for j = i+1:N,
%                 phi2 = phi(j);
%                 if A(i,j) ~= 0                    
%                     polarplot([phi1 phi2],[1 1],':b')
%                     hold on
%                     avg = (phi1 + phi2)/2;
%                     if abs(phi2 - phi1) > pi
%                         avg = avg + pi;
%                         ravg = cos(avg-phi2);
%                     else
%                         ravg = cos(avg-phi1);
%                     end
%                     text(avg,1.1*ravg,num2str(A(i,j)))
%                 end
%             end
%             if i < N,
%                 polarplot([phi(B(i,2)) phi(B(i,3))],[1 1],'r')
%             end
%         end
        % Nodes
        polarscatter(phi,rho,400,'filled','MarkerEdgeColor','k','MarkerFaceColor','w')
        for i = 1:N,
            text(phi(i),1,num2str(i),'HorizontalAlignment','Center');
        end
        set(gca,'ThetaTickLabel',[])
        set(gca,'RTickLabel',[])
        set(gca,'RLim',[0 1.2])
        grid off
        set(gcf,'Color','w')
        title('Spanning Tree Construction')
        dim = [0.8 0.8 0.1 0.1];
        annotation('textbox',dim,'String',{'Original',''},'Color','b','BackgroundColor','w')
        annotation('textbox',dim,'String',{'','MST'},'Color','r','LineStyle','none')
    else
        varargout{1} = B;
        if nargout == 2,
            varargout{2} = wt;
        end
    end
    
    
    % New find function
    function y = sfind(parent,idx)
        if parent(idx) ~= idx,
            y = sfind(parent, parent(idx));
        else
            y = idx;
        end
    end
    
    % New union function
    function [pout,rout] = subunion(parent,rank,i,j)
        xrt = sfind(parent,i);
        yrt = sfind(parent,j);
        
        if rank(xrt) < rank(yrt)
            parent(xrt) = yrt;
        elseif rank(xrt) > rank(yrt)
            parent(yrt) = xrt;
        else
            parent(yrt) = xrt;
            rank(xrt) = rank(xrt) + 1;
        end
        pout = parent;
        rout = rank;
    end
    
end

