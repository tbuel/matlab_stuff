function [rfV, rfH] = RF_process(frfvert,frfhorz,option)
% function [rfV, rfH] = RF_process(frfvert,frfhorz,option)
%   RF_process takes the RF pattern images and retrieves data for use in 
%   plotting and other calculations. RF_process also plots a comparison
%   between the input images, the reconstructed polar plots, and the
%   reconstructed cartesian plots.
%
%   - frfvert = filename for RF vertical pattern image file
%   - frfhorz = filename for RF horizontal pattern image file
%   - rfV = two-column matrix, of the vertical RF pattern radiation. The 
%           1st column is the angle in radians, and the 2nd column is the 
%           normalized unitless gain. The gain has values between 0 and 1, 
%           where 1 = 0dB.
%   - rfV = two-column matrix, of the horizontal RF pattern radiation. The 
%           1st column is the angle in radians, and the 2nd column is the 
%           normalized unitless gain. The gain has values between 0 and 1, 
%           where 1 = 0dB.
%
%   Additionally, if no filenames are given, the program will prompt the
%   user for the filenames. The program also prompts the user for radial
%   datum points on the original RF pattern images. If these files already
%   exist (As {frfvert}_db.txt and {frfhorz}_db.txt), the program uploads them.
%   Be aware of naming any text files with that description!
%
%   Tanner Buel 12/22/2017
%
%% ------------ System Characteristics ------------
c = 299792458; % Speed of light; c = ~3*10^8 m/s

%% ------------ Read RF Antenna Patterns ------------
theta = 0:0.01:2*pi;
if nargin == 0,
    fprintf('RF Process\n');
    frfvert = input('Enter RF vertical image file: ','s');
    if isempty(frfvert),
        frfvert = 'Ref\HG5158P-RF_vert.JPG';
    end
    frfhorz = input('Enter RF horizontal image file: ','s');
    if isempty(frfhorz),
        frfhorz = 'Ref\HG5158P-RF_horz.JPG';
    end
end
frfvshort = frfvert(1:length(frfvert)-4);
frfhshort = frfhorz(1:length(frfhorz)-4);
Vpol = imread(frfvert);
Hpol = imread(frfhorz);

% pick red points
Vpts = Vpol(:,:,1) > 127 & Vpol(:,:,2) < 127 & Vpol(:,:,3) < 127; 
Hpts = Hpol(:,:,1) > 127 & Hpol(:,:,2) < 127 & Hpol(:,:,3) < 127;
[xv, yv] = find(Vpts);
[xh, yh] = find(Hpts);
% [xv, yv] = meshgrid(1:size(Vpts,1), 1:size(Vpts,2));
% [xh, yh] = meshgrid(1:size(Hpts,1), 1:size(Hpts,2));

% get scale
% % vopt = input('Enter filename for vertical RF db points:','s');
out3db = zeros(6,2);
out3db(:,1) = [0,3,10,20,30,40]';
figure(1)
image(Vpol)
title('Vertical RF Plot')
fprintf('\nClick on center of vertical RF plot\n');
[xvc,yvc] = ginput(1);
fRvdb = [frfvshort '_db.txt'];
if ~exist(fRvdb,'file')
    fprintf('\nClick on circles at R = [0db, -3db, -10db, -20db, -30db, -40db] in order, on vertical RF plot\n');
    [xsv,ysv] = ginput(6);
    Rvdb = sqrt((xsv-xvc).^2 + (ysv-yvc).^2);
    Rvdb = (Rvdb./max(Rvdb));%.*ones(length(theta),1);
    out3db(:,2) = Rvdb;
    save(fRvdb,'-ascii','-tabs','out3db');
else
    Rvdb = load(fRvdb);
    Rvdb = Rvdb(:,2);
end
cla
image(Hpol)
title('Horizontal RF Plot')
fprintf('\nClick on center of horizontal RF plot\n');
[xhc,yhc] = ginput(1);
fRhdb = [frfhshort '_db.txt'];
if ~exist(fRhdb,'file')
    fprintf('\nClick on circles at R = [0db, -3db, -10db, -20db, -30db, -40db] in order, on horizontal RF plot\n');
    [xsh,ysh] = ginput(6);
    Rhdb = sqrt((xsh-xhc).^2 + (ysh-yhc).^2);
    Rhdb = (Rhdb./max(Rhdb));%.*ones(length(theta),1);
    out3db(:,2) = Rhdb;
    save(fRhdb,'-ascii','-tabs','out3db');
else
    Rhdb = load(fRhdb);
    Rhdb = Rhdb(:,2);
end
close(1)




% ------------ Scale and Normalize RF Pattern ------------
% xh = xh - mean(xv);
% yv = yv - size(Vpts,2)/2;
% xh = xh - mean(xh);
% yh = yh - size(Hpts,2)/2;

xv = xv - mean(xv);
yv = yv - yvc; %size(Vpts,2)/2;
xh = xh - mean(xh);
yh = yh - yhc; %size(Hpts,2)/2;


% Convert to Polar coordinates
[phiv,rv] = cart2pol(xv,yv);
[phih,rh] = cart2pol(xh,yh);
rv = rv./max(rv);
rh = rh./max(rh);
phiv = phiv - pi/2;
phih = phih - pi/2;


%% ------------ Plot RF Antenna Patterns ------------
% Original Images
figure
set(gcf,'color','w')
subplot(2,3,1)
imshow(Vpol)
title('Vertical RF Pattern - Imported')
subplot(2,3,4)
imshow(Hpol)
title('Horizontal RF Pattern - Imported')

% Reconstructed
phitickval = 0:10:360;
rvtickval = flip(Rvdb);%[1/16 1/8 1/4 1/2 9/10 1];
rticklab = {'-40','-30','-20','-10','-3','0'};
subplot(2,3,2)
angle = 1;
while ~isempty(angle),
	polarplot(phiv,rv,'.','MarkerSize',1);
    title('Vertical RF Pattern - Recontructed');
	angle = input('\nEnter adjustment angle [deg] for vertical pattern: ');
    if ~isempty(angle)
        phiv = phiv + angle*pi/180;
    end
    if angle == 0,
        break;
    end
end
hold on
axv = gca;
axv.ThetaMinorGrid = 'on';
axv.GridAlpha = 0.4;
axv.ThetaAxis.MinorTickValues = phitickval;
axv.ThetaAxis.LineWidth = 2;
axv.RTick = rvtickval;
axv.RTickLabel = rticklab;
axv.RAxisLocation = 45;

subplot(2,3,5)
rhtickval = flip(Rhdb);
angle = 1;
while ~isempty(angle),
	polarplot(phih,rh,'.','MarkerSize',1);   
    title('Horizontal RF Pattern - Recontructed');
	angle = input('\nEnter adjustment angle [deg] for horizontal pattern: ');
    if ~isempty(angle)
        phih = phih + angle*pi/180;
    end
    if angle == 0,
        break;
    end

end
axh = gca;
axh.ThetaMinorGrid = 'on';
axh.GridAlpha = 0.4;
axh.ThetaAxis.MinorTickValues = phitickval;
axh.ThetaAxis.LineWidth = 2;
axh.RTick = rhtickval;
axh.RTickLabel = rticklab;
axh.RAxisLocation = 45;

% Filter, smooth, shifting
if length(phiv) == 1,
    phiv = phiv';
    rv = rv';
    phih = phih';
    rh = rh';
end    
for p = 1:length(phiv),
    if phiv(p) < -pi,
        phiv(p) = phiv(p) + 2*pi;
    end
end
for p = 1:length(phih),
    if phih(p) < -pi,
        phih(p) = phih(p) + 2*pi;
    end
end
rfV = [phiv rv];
rfH = [phih rh];
rfV = sortrows(rfV,1);
rfH = sortrows(rfH,1);
save([frfvshort '.txt'],'-ascii','-tabs','rfV');
save([frfhshort '.txt'],'-ascii','-tabs','rfH');
phiv = rfV(:,1);
rv = rfV(:,2);
phih = rfH(:,1);
rh = rfH(:,2);
rv = movavg(rv,10);
rh = movavg(rh,10);


%figure
subplot(2,3,3)
plot(phiv*180/pi,rv)
title('Vertical RF Pattern - Recontructed - Cartesian')
xlabel('Angle, deg')
ylabel('Gain, dbi')
xlim([-180,180])
xticks(-180:30:180)
yticks(rvtickval)
yticklabels(rticklab)
grid on

subplot(2,3,6)
plot(phih*180/pi,rh)
title('Horizontal RF Pattern - Recontructed - Cartesian')
xlabel('Angle, deg')
ylabel('Gain, dbi')
xlim([-180,180])
xticks(-180:30:180)
yticks(rhtickval)
yticklabels(rticklab)
grid on

%% ----------------- EMR Plots --------------------------
if nargin > 2
    if option == 'show3d'
        [xv,yv] = pol2cart(phiv,rv);
        [xh,yh] = pol2cart(phih,rh);
        figure
        scatter3(xv,zeros(size(xv)),yv,'.')
        hold on
        scatter3(xh,yh,zeros(size(xh)),'.')
        box on

        xc = 0.5;
        yc = 0;
        zc = median(yv);
        xR = max(xv)/2;
        yR = 0.58;
        zR = 0.5;
        Ne = 20;
        [xe, ye, ze] = ellipsoid(xc,yc,zc,xR,yR,zR,Ne);
        surf(xe,ye,ze)
        title('EM Radiation Approximation')
        xlabel('x-axis')
        ylabel('y-axis')
        zlabel('z-axis')
        legend('vpol','hpol')
    end
end
end
