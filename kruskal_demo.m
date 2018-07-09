%% Kruskal Algorithm Demo
% TB 1/11/18
close all
clearvars
clc

% User Input Spanning Tree Network
opt = input('\nEnter number of vertices of graph (or enter to generate random graph): ');
if ~isempty(opt) 
    N = min(opt,10);
    A = zeros(N);
    for i = 1:N,
        row = input('\nEnter %dth row: ');
        A(i,:) = row;
    end
else
    % Generate Random Network (Very Busy)
    N = randi([4 10],1); % number of nodes/vertices
    A = zeros(N);
    lim = 20;
    for i = 1:N,
        x = randi([0 lim],N-i,1);
        A(i,i+1:N) = x';
    end
    A = A + A';
end
fprintf('\nNumber of Vertices: %d\n',N);
fprintf('\nOriginal Tree:\n');
disp(A);
wt0 = sum(sum(A))/2; % original weight

% Minimum Spanning Tree - Kruskal's Algorithm
[B,wt] = kruskal(A);

% Plotting - This code is generated if kruskal(A) was not assigned an
% output
fprintf('\nMinimum Spanning Tree:\n\tWt\tNode1\tNode2\n');
disp(B);
fprintf('Input Tree Weight: %d\nMST Weight: %d\nSaved Weight: %.2f%%\n',wt0,wt,100*(wt0-wt)/wt0);

figure
phi = (0:N-1)*2*pi/N;
rho = ones(1,N);
% Spanning Tree
for i = 1:N,
    phi1 = phi(i);
    for j = i+1:N,
        phi2 = phi(j);
        if A(i,j) ~= 0                    
            polarplot([phi1 phi2],[1 1],':b')
            hold on
            avg = (phi1 + phi2)/2;
            if abs(phi2 - phi1) > pi
                avg = avg + pi;
                ravg = cos(avg-phi2);
            else
                ravg = cos(avg-phi1);
            end
            text(avg,1.1*ravg,num2str(A(i,j)))
        end
    end
    if i < N,
        polarplot([phi(B(i,2)) phi(B(i,3))],[1 1],'r')
    end
end
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