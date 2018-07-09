function mass_spring_ode(m,k,c)
% function mass_spring_ode(mass,kstiff,nudamp)
%   Testing out Matlab's ODE functions with a mass on a spring

close all
if nargin == 0,
    m = 5; %kg
    k = 100; 
    c = 5;
end
tspan = [0 10];
g = -9.81; %acceleration due to gravity
y0 = -5;
dy0 = 0;
yf = 0;
[t,y] = ode45(@f ,tspan, [y0 dy0]);
dt = t(2)-t(1);


figure(1)
hx1 = subplot(1,2,1);
ax1 = plot(t,y(:,1),'-o');
pBoxR = get(hx1,'plotboxaspectratio')
outPos = get(hx1,'outerposition')
outPos(3) = outPos(3) + 0.1;
outPos(1) = outPos(1) - 0.05;
set(hx1,'outerposition',outPos)
plotpos1 = get(hx1,'position');
plotpos1(3) = plotpos1(3)+0.05;
set(hx1,'position',plotpos1)
xlabel('Time, s')
ylabel('Height, m')
title('Mass Spring Damper System')
ysym(hx1);
ylim1 = ylim;
grid on

N = length(t);
% figure
hx2 = subplot(1,2,2);
plotpos2 = get(hx2,'position');
plotpos2(3) = 0.05;
plotpos2(1) = plotpos1(1) + plotpos1(3) + 0.2;
set(hx2,'position',plotpos2);
for i=1:N,
    ax2 = stem(hx2,0,y(i,1));
    xlim([-1 1])
    set(hx2,'XTickLabel',[])
    set(hx2,'YTickLabel',[])
    title('Animation')
    ylim(ylim1)
    pause(dt)
end


    function dydt = f(t,y)
        dydt = [y(2); - (k/m)*y(1) - (c/m)*y(2) + g];
    end
end

