%% Relative Ball Drop
function reldrop(Ht,parray,option)
    
    if nargin == 0,
        Ht = 5;
        parray = [];
        option = 'init';
    elseif nargin <3,
        option = 'init';
    end
    if isempty(parray)
        parray = 1:9;
    end
    splanets = {'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};
    rGrav = 4; % Gravity
    fps = 15; % frames per second for "videos"
    ts = 1/fps;
    if ~exist('data','var')
        data = dlmread('planets.txt','\t',1,1);
    end
    
    if strcmp(option,'init'),
        hf = figure('Name','Relative Ball Drop on each Planet','NumberTitle','off');
    elseif strcmp(option,'replot'),
        hf = gcf;
        cla;
    end
    splotwid = 0.4;
    sploty = 0.4;
    splotht = 0.5;
    hx1 = subplot(1,2,1)
    set(hx1,'position',[0.1 sploty splotwid splotht])
    hold on
    steps = 2^6; % Number of steps for 1 Revolution
    ax1 = gca;
    hts = Ht*ones(1,9);
    dt = ts; %delta time interval
    addlabels = 0;
    vlast = zeros(1,9); % initial condition velocity m/s
    vy = zeros(1,9); % initial condition velocity m/s
    vf = zeros(1,9); % final velocity
    vexp = sqrt(2.*data(rGrav,:)*Ht);
    t = 0;
    cols = [];
    ii = 1;
    while ~isequal(hts,zeros(1,9))
        if exist('h2','var')
            delete(h2)
        end
        h2 = gobjects(1,9);
        for n = 1:9,
            if parray(n) ~= n,
                continue;
            else
                cols(ii) = n;
                ii = ii + 1;
            end
            h2(n) = plot(ax1,n-1,hts(n),'.b','MarkerSize',20);
            vy(n) = vlast(n) - data(rGrav,n)*dt;
            hts(n) = hts(n) + vy(n)*dt;
            vlast(n) = vy(n);
            if hts(n) <= 0
                hts(n) = 0;
                if vf(n) == 0,
                    vf(n) = vlast(n);
                end
            end
            axis([-1 9 0 Ht])
            if addlabels == 0,
                ylabel('Height, m')
                xlabel('Planet')
                addlabels = 1;
                grid on
                xticks(0:8);
                xticklabels(splanets);
                xtickangle(60);
            end
        end
        telapse = ['Time Elapsed: ' num2str(t,'%.2f') 's'];
        stitle = {'Drop Comparison on each Planet'; telapse};
        title(stitle)
        t = t + dt;
        pause(ts/10)
    end
    tf = t;

    t = 0:dt:tf;
    N = length(t);
    xy = Ht - data(rGrav,:).*t'.^2/2;
    hx2 = subplot(1,2,2)
    ax2 = gca;
    set(hx2,'position',[0.15 + splotwid sploty splotwid splotht])
    plot(ax2,t,xy(:,cols))
    title('Time History')
    ylabel('Height, m')
    xlabel('Time, s')
    ylim([0 Ht])
    grid on
    legend(splanets)
    
    
    %% UI Control
    % Choose which planets to see
    chkbx = zeros(9,1);
    for ck = 1:9,
        chkbx(ck) = uicontrol(hf,'Style', 'checkbox',...
                           'String', splanets{ck},...
                           'Units','normalized', ...
                           'Position',  [(ck*0.1 - 0.05) 0.2 0.15 .05],...
                           'Value',1); 
    end
   set(chkbx,'Callback',{@setplanet,chkbx});       
       
    %Show Again
    btn1 = uicontrol(hf,'Style', 'pushbutton',...
           'String','Show!',...
           'Units','normalized', ...
           'Position',  [0.45 0.05 0.1 0.1],...
           'Callback',@reshow); 
       
    % Checkbox Callback
    function setplanet(hObject,eventData,checkBoxID)
        value = get(hObject,'Value');
        if value 
            fprintf('\n %s',splanets(checkBoxID));
            switch checkBoxID
                case 1
                    fprintf('');
                case 2
                    fprintf('');
                case 3
                    fprintf('');
                case 4
                    fprintf('');
                case 5
                    fprintf('');
                case 6
                    fprintf('');
                case 7
                    fprintf('');
                case 8
                    fprintf('');
                case 9
                    fprintf('');
                otherwise
                    fprintf('Invalid\n');
            end
        end
%             newPlanet = get(popup1,'Value')
    end

    % Slider Callback
    function setrev(src,event)
         newRev = get(slider1,'Value');
         set(txt3, 'String', num2str(newRev));
    end

    % Button Callback
    function reshow(src,event)
         reldrop(get(popup1,'Value'),get(slider1,'Value'),'replot');
    end

end