N = 2^6;
phi = linspace(0,2*pi,N);
x = cos(phi);
y = sin(phi);
p = 0;
dt = 1/30; % pause interval
r = 10; % rate
dp = phi(2) - phi(1);
while(1)
    p = p + dp;
    g = 1;%round(exp(-p/r),2);
    z = g*cos(phi + p);
    z = z - min(z);
    H = stem3(x,y,z,'^');
%     H = plot3(x,y,z,'-o');
    grid on
    zlim([0 2]);
    pause(dt);
    if max(z) < 0.001 || ~ishandle(H)
        break;
    end
end