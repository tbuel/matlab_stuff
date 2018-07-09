function orbit(aPlanet,NRev,option)
% Show Orbits of planets
% aPlanet must be an integer 1-9 (yea I know pluto isn't a planet)
% NRev must be an integer 1-5
% Option to come
% tbuel 1/19/18
    if nargin < 2
        aPlanet = 3;
        NRev = 1;
    end
    if nargin < 3
        option = 'init';
    end
    
    G = 6.6741e-11; % Nm^2/kg^2
    N = 2^12;
    phi = (0:N-1)*(2*pi)/N;
    fps = 15; % frames per second for "videos"
    ts = 1/fps;
    splanets = {'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};
    rOrb = 11;
    if ~exist('data','var')
        data = dlmread('planets.txt','\t',1,1);
    end


    %% Orbits - Assume circular orbit for now
    pos = [1 0 0; 2 0 0; 3 0 0; 4 0 0; 5 0 0; 6 0 0; 7 0 0; 8 0 0; 9 0 0];
    plotpos = [0.05 0.1 0.6 0.8]; 
    if strcmp(option,'init'),
        hf = figure('Name','Orbits',...
            'NumberTitle','off');
    elseif strcmp(option,'replot'),
        hf = gcf;
        cla;
    end
    hx = subplot(1,2,1); % Use Subplot for annotation at right
    plot(0,0,'.y','MarkerSize',50) % Plot Sun
    title('Solar System')
    ax2 = gca;
    set(ax2,'color','k')
    set(hx,'position',plotpos)
    hold on
    % Orbital Paths
    for n = 1:9,
        xOrb = n*cos(phi);
        yOrb = n*sin(phi);
        plot(ax2,xOrb,yOrb,'-w')            
    end
    set(ax2,'XTickLabel',[]);
    set(ax2,'YTickLabel',[]);
    sub2x = (plotpos(1) + plotpos(3) + 0.05); % position of information

    pcols = [[128/255,128/255,128/255];[204/255, 102/255, 0];[0 0 1];[1 0 0];[1, 163/255, 26/255];[194/255, 151/255, 10/255];[102/255, 1, 204/255];[0, 51/255, 204/255];[0 1 1]];
    NperRev = 2^6; % Number of steps for 1 Revolution
    Tyr = data(rOrb,aPlanet); % 1 Revolution at Relative Planet
    inc = Tyr*2*pi/NperRev; % degree increment
    boxdim = [sub2x (plotpos(2) + plotpos(4) - 0.1) .1 .1]; % annotation dimensions
    for rev = 1:NRev,
        for i = 1:NperRev,
            if exist('harray','var')
                delete(harray)
            end
            harray = gobjects(1,9);
            for n = 1:9,
                harray(n) = plot(ax2,pos(n,1),pos(n,2),'.','Color',pcols(n,:),'MarkerSize',20);
                theta = atan2(pos(n,2),pos(n,1));
                r = sqrt(pos(n,1)^2 + pos(n,2)^2);
                theta = mod(theta + inc/data(rOrb,n), 2*pi);
                pos(n,1) = r*cos(theta);
                pos(n,2) = r*sin(theta);
                pos(n,3) = (pos(n,3) + Tyr/NperRev/data(rOrb,n));
            end
            planetyears = { 'Number of Revolutions:';
                            ['Mercury: ' num2str(pos(1,3),'%.4g')];...
                            ['Venus: ' num2str(pos(2,3),'%.4g')];...
                            ['Earth: ' num2str(pos(3,3),'%.4g')];...
                            ['Mars: ' num2str(pos(4,3),'%.4g')];...
                            ['Jupiter: ' num2str(pos(5,3),'%.4g')];...
                            ['Saturn: ' num2str(pos(6,3),'%.4g')];...
                            ['Uranus: ' num2str(pos(7,3),'%.4g')];...
                            ['Neptune: ' num2str(pos(8,3),'%.4g')];...
                            ['Pluto: ' num2str(pos(9,3),'%.4g')]};
            if exist('ta','var'),
                delete(ta);
            end
            ta = annotation(hf,'textbox',boxdim,'String',planetyears,'FitBoxToText','on','BackgroundColor','w');
            pause(ts)
        end
    end
    hold off
    
    %% UI Control
    uiyoff = ta.Position(2) - 0.1;
    uiwid = ta.Position(3);
    % Get Planets
    txt1 = uicontrol(hf,'Style', 'text',...
           'String', 'Select Planet:',...
           'Units','normalized', ...
           'HorizontalAlignment','Left',...
           'Position',  [sub2x uiyoff uiwid/2 .05]);  
    popup1 = uicontrol(hf,'Style', 'popup',...
           'String', splanets,...
           'Units','normalized', ...
           'Position',  [sub2x+uiwid/2 uiyoff uiwid/2 .05],...
           'Value',aPlanet,...
           'Callback',@setplanet);   
    % Get Revolutions
    txt2 = uicontrol(hf,'Style', 'text',...
           'String', 'New # of Revolutions:',...
           'Units','normalized', ...
           'HorizontalAlignment','Left',...
           'Position',  [sub2x (uiyoff - 0.1) 0.9 .05]);  
    txt3 = uicontrol(hf,'Style', 'text',...
           'String', num2str(NRev),...
           'Units','normalized', ...
           'HorizontalAlignment','Left',...
           'Position',  [(sub2x+0.9*uiwid) (uiyoff - 0.1) 0.1 .05]);  
    slider1 = uicontrol(hf,'Style', 'slider',...
           'Min',1,'Max',5,'Value',NRev,...
           'SliderStep',[0.25 0.25],...
           'Units','normalized', ...
           'Position',  [sub2x (uiyoff - 0.15) uiwid .05],...
           'Callback',@setrev); 
    %Show Again
    btn1 = uicontrol(hf,'Style', 'pushbutton',...
           'String','Show!',...
           'Units','normalized', ...
           'Position',  [(sub2x + uiwid/2 - 0.05) 0.1 0.1 0.1],...
           'Callback',@reshow); 
    
    % Drop-down menu Callback
    function setplanet(src,event)
%         newPlanet = get(popup1,'Value')
    end

    % Slider Callback
    function setrev(src,event)
         newRev = get(slider1,'Value');
         set(txt3, 'String', num2str(newRev));
    end

    % Button Callback
    function reshow(src,event)
         orbit(get(popup1,'Value'),get(slider1,'Value'),'replot');
    end
end

