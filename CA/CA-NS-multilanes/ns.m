function [rho, flux, vmean] = ns(rho, p, L, tmax, animation, spacetime)
%
% NS: This script implements the Nagel Schreckenberg cellular automata based
% traffic model. Car move forward governed by NS algorithm:
% 
%   1. Acceleration. If the vehicle can speed up without hitting the speed 
%      limit vmax it will add one to its velocity, vn -> vn + 1. Otherwise, 
%      the vehicle has constant speed, vn -> vn.
%
%   2. Collision prevention. If the distance between the vehicle and the car
%      ahead of it, dn, is less than or equal to vn , i.e. the nth vehicle
%      will collide if it doesn't slow down, then vn -> dn − 1.
%
%   3. Random slowing. Vehicles often slow for non-traffic reasons (cell 
%      phones, coffee mugs, even laptops) and drivers occasionally make 
%      irrational choices. With some probability pbrake, vn -> vn − 1, 
%      presuming vn > 0.
%
%   4. Vehicle movement. The vehicles are deterministically moved by their 
%      velocities, xn -> xn + vn .
%
% USAGE: flux = ns(rho, p, L, tmax, isdraw)
%        rho       = density of the traffic
%        p         = probability of random braking
%        L         = length of the load
%        tmax      = number of the iterations
%        animation = if show the animation of the traffic
%        spacetime = if plot the space-time after the simuation ended.
%        flux      = flux of the traffic
%
% zhou lvwen: zhou.lv.wen@gmail.com
%

if nargin == 0; 
    rho = 0.15; p = 0.25; L = 100; tmax = 5000; pchange = 0.5;
    animation = 'circle';
end

%rand('seed',1)
vmax = 5;                      % maximun speed
% place a distribution with density
ncar = round(L*2*rho);
rho = ncar/2/L;

xy = randperm(2*L,ncar);
[y,x] = ind2sub([2, L], xy); % y: 2 for right, 1 for left

v = vmax * ones(1,ncar);       % start everyone initially at vmax

switch animation
    case 'circle'; h = plotcirc(L,x,y,2);
    case 'line';   h = plotline(L,x,y,2);
end

flux = [0 0];                  % number of cars that pass through the end
vmean = [0 0];

voffset = 1;
vback = 1;
XY = [x y];
for t = 1:tmax
    % determine the space vehicles have to move
    [gaps, gapfront, gapback] = gaplength(x,y,L); 
    
    % left to right & right to left
    l2r = find(y==1 & gaps>vmax+voffset & gapfront>vmax+voffset & gapback>=vback);
    r2l = find(y==2 & gaps<vmax         & gapfront>gaps         & gapback>=vback & rand(size(y))<pchange);
    y(l2r) = 2;
    y(r2l) = 1;
    % acceleration
    v = min(v+1, vmax);
    
    gaps = gaplength(x,y,L); 
    %collision prevention
    v = min(v, gaps-1);
    
    % random speed drops
    vdrops = ( rand(1,ncar)<p );
    v = max(v-vdrops,0);
    
    % update the position
    x = x + v;
    passed = x>L;             % cars passed at time r
    x(passed) = x(passed) - L;% periodic boundary conditions
    
    if t>tmax/2
        flux(1) = flux(1) + sum(v(y==1)/L); %flux = flux + sum(passed); 
        flux(2) = flux(2) + sum(v(y==2)/L); %flux = flux + sum(passed); 
        vmean(1) = vmean(1) + mean(v(y==1));
        vmean(2) = vmean(2) + mean(v(y==2));
    end
    
    switch animation
        case 'circle'; plotcirc(L,x,y,0.1,h); 
        case 'line';   plotline(L,x,y,0.1,h);
    end
    
end
flux = flux/(tmax/2);
vmean = vmean/(tmax/2);

% -------------------------------------------------------------------------


function [gap, gapfront, gapback] = gaplength(x,y,L)
% 
% GAPLENGTH: determine the gaps between vehicles
%
ncar = length(x);
gap = inf*ones(1, ncar);
gapfront = inf*ones(1, ncar);
gapback = inf*ones(1, ncar);
index = 1:ncar;
for i = index
    j1 = index(index~=i & y==y(i));
    if ~isempty(j1)
        d1 = x(j1) - x(i);
        d1(d1<-L/2) = d1(d1<-L/2) + L;
        if any(d1>0)
            gap(i) = min(d1(d1>0));
        end
    end
    
    j2 = index(index~=i & y~=y(i));
    if ~isempty(j2)
        d2 = x(j2) - x(i);
        d2(d2<-L/2) = d2(d2<-L/2) + L;
        if any(d2>=0)
            gapfront(i) = min(d2(d2>=0));
        end
        d3 = x(i) - x(j2);
        d3(d3<-L/2) = d3(d3<-L/2) + L;
        if any(d3>=0)
            gapback(i)  = min(d3(d3>=0));
        end
    end
end

% -------------------------------------------------------------------------

function h = plotcirc(L,x,y,dt,h)
W = 0.05;  
ncar = length(x);

theta = [(0-pi/L) : 2*pi/L : (2*pi+pi/L)];
R = ones(size(theta));
theta = [ theta; theta];
R = [R; R+W];
xc = cos(theta);     yc = sin(theta);
xinner = (R-W/2).*xc; yinner = (R-W/2).*yc;
xouter = (R+W/2).*xc; youter = (R+W/2).*yc;
i = sub2ind(size(R),y,x);
if nargin == 4
    color = randperm(ncar);
    xi = [xinner(i);  xinner(i+2); xouter(i+2); xouter(i)];
    yi = [yinner(i);  yinner(i+2); youter(i+2); youter(i)];
    h = fill(xi,yi, color); hold on
    plot(xinner(1,:),yinner(1,:), 'k', ...
         xouter(1,:),youter(1,:), 'k', ...
         xouter(2,:),youter(2,:), 'k','linewidth',1.5);
    plot([xinner; xouter], [yinner; youter],'k','linewidth',1.5);
    axis image; 
else
    xi = [xinner(i);  xinner(i+2); xouter(i+2); xouter(i)];
    yi = [yinner(i);  yinner(i+2); youter(i+2); youter(i)];
    for i=1:ncar;  set(h(i),'xdata',xi(:,i),'ydata',yi(:,i));  end
end
pause(dt)

% -------------------------------------------------------------------------

function h = plotline(L,x,y,dt,h)
W = 2;
rmin = 1;
dw = 0.05; 
rmax = rmin + dw*(W-1);
ncar = length(x);
ti = 0 : 2*pi/L : 2*pi;
ri = rmin:dw:rmax;

[theta, R] = meshgrid(ti, ri);

xmin = (R-dw/2).*cos(theta);
ymin = (R-dw/2).*sin(theta);
xmax = (R+dw/2).*cos(theta);
ymax = (R+dw/2).*sin(theta);

i = sub2ind(size(R),y,x);
xi = [xmin(i);  xmin(i+W); xmax(i+W); xmax(i)];
yi = [ymin(i);  ymin(i+W); ymax(i+W); ymax(i)];

if nargin == 6
    color = randperm(ncar);
    h = fill(xi,yi, color); hold on
    plot([xmin; xmax]', [ymin; ymax]', 'k', 'linewidth', 1.5);
    plot([xmin; xmax] , [ymin; ymax] , 'k', 'linewidth', 1.5);
    axis image;
else
    for i=1:ncar;  set(h(i),'xdata',xi(:,i),'ydata',yi(:,i));  end
end
pause(dt)
