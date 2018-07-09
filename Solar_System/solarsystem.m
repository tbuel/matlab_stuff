%% Solar System Stuff
% TB 1/4/2018
clearvars -except data
close all
clc

%% Constants
G = 6.6741e-11; % Nm^2/kg^2
N = 2^12;
phi = (0:N-1)*(2*pi)/N;
fps = 30; % frames per second for "videos"
ts = 1/fps;
splanets = {'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};

% Load from file
% Column Numbers
rMass = 1; % Mass
rDia = 2; % Diameter
rDens = 3; % Density
rGrav = 4; % Gravity
rEV = 5; % Escape Velocity
rRP = 6; % Rotation Period
rDay = 7; % Length of Day
rDist = 8; % Distance from Sun
rPeri = 9; % Perihelion
rAph = 10; % Aphelion
rOrb = 11; % Orbital Period (days)
rOV = 12; % Orbital Velocity (km/s)
rOI = 13; % Orbital Inclination
rOE = 14; % Orbital Eccentricity
rObl = 15; % Obliquity to Orbit (deg)
rTemp = 16; % Temp, C
rMoon = 17; % N moons
rRing = 18; % Rings?
rHField = 19; % Magnetic Field
data = dlmread('planets.txt','\t',1,1);
data(rMass,:) = data(rMass,:)*1e24; % Mass (10^24kg)
data(rDist,:) = data(rDist,:)*1e6; % Distance from Sun (10^6km)
data(rPeri,:) = data(rPeri,:)*1e6; % Perihelion (10^6km)
data(rAph,:) = data(rAph,:)*1e6; % Aphelion(10^6km)


% Planets
Mercury = data(:,1);
Venus = data(:,2);
Earth = data(:,3);
Mars = data(:,4);
Jupiter = data(:,5);
Saturn = data(:,6);
Uranus = data(:,7);
Neptune = data(:,8);
Pluto = data(:,9);


%% Scale Image
f1 = figure;
plot(0,0,'.y','MarkerSize',2*data(rDia,5)/data(rDia,9))
title('Solar System to Scale - Except Sun')
ax = gca;
set(ax,'color','w')
hold on
for n = 1:9,
    R = data(rDist,n);
    dia = data(rDia,n)/data(rDia,9);
    plot(R,0,'.','MarkerSize',dia)
    a = data(rAph,n)/(1+data(rOE,n));
    b = a*sqrt(1-data(rOE,n).^2);
    xOrb = a*cos(phi)+data(rOE,n);
    yOrb = b*sin(phi);
    plot(xOrb,yOrb,'-k')
end
hold off

%% Orbits
% Get Planet
aPlanet = 3;
while ~(aPlanet >= 1 && aPlanet <= 9)
    aPlanet = input('\nSelect planet, [1-9], to see relative orbital revolutions: ');
end

% Get Number of REvolutions
NRev = 1;
while ~(NRev >= 1 && NRev <= 5)
    NRev = input('\nEnter number of revolutions [1-5]: ');
end
orbit(aPlanet,NRev);


%% Relative Ball Drop
aPlanet = 3;
while ~(aPlanet >= 1 && aPlanet <= 9)
    aPlanet = input('\nSelect planet, [1-9], to see relative ball drop: ');
end
Ht = 5;
while ~(Ht >= 1 && Ht <= 10)
    Ht = input('\nEnter height from ground: ');
end
reldrop(Ht,[]);
