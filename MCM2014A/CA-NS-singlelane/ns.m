function [rho, flux, vmean] = ns(rho, p, L, tmax, animation, spacetime)
%
% NS: This script implements the Nagel Schreckenberg cellular automata based
% traffic model. Car move forward governed by NS algorithm:
% 
%   1. Acceleration. If the vehicle can speed up without hitting the speed 
%      limit vmax it will add one to its velocity, vn → vn + 1. Otherwise, 
%      the vehicle has constant speed, vn → vn.
%
%   2. Collision prevention. If the distance between the vehicle and the car
%      ahead of it, dn, is less than or equal to vn , i.e. the nth vehicle
%      will collide if it doesn't slow down, then vn → dn − 1.
%
%   3. Random slowing. Vehicles often slow for non-traffic reasons (cell 
%      phones, coffee mugs, even laptops) and drivers occasionally make 
%      irrational choices. With some probability pbrake, vn → vn − 1, 
%      presuming vn > 0.
%
%   4. Vehicle movement. The vehicles are deterministically moved by their 
%      velocities, xn → xn + vn .
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
    rho = 0.25; p = 0.25; L = 100; tmax = 100; 
    animation = true; spacetime = true;
end

vmax = 5;                  % maximun speed
% place a distribution with density
ncar = round(L*rho);
rho = ncar/L;

x = sort(randsample(1:L, ncar));
v = vmax * ones(1,ncar);   % start everyone initially at vmax

if animation; h = plotcirc(L,x,0.1); end

flux = 0;                  % number of cars that pass through the end
vmean = 0;
road = zeros(tmax, L);

for t = 1:tmax
    % acceleration
    v = min(v+1, vmax);
    
    %collision prevention
    gaps = gaplength(x,L); % determine the space vehicles have to move
    v = min(v, gaps-1);
    
    % random speed drops
    vdrops = ( rand(1,ncar)<p );
    v = max(v-vdrops,0);
    
    % update the position
    x = x + v;
    passed = x>L;          % cars passed at time r
    x(passed) = x(passed) - L;% periodic boundary conditions
    
    if t>tmax/2
        flux = flux + sum(v/L); %flux = flux + sum(passed); 
        vmean = vmean + mean(v);
    end
    road(t,x) = 1;
    
    if animation; h = plotcirc(L,x,0.1,h); end
end
flux = flux/(tmax/2);
vmean = vmean/(tmax/2);

if spacetime; figure;imagesc(road);colormap([1,1,1;0,0,0]);axis image; end

% -------------------------------------------------------------------------

function gaps = gaplength(x,L)
% 
% GAPLENGTH: determine the gaps between vehicles
%
ncar = length(x);
gaps=zeros(1, ncar);
if ncar>0
    gaps = x([2:end 1]) -x;
    gaps(gaps<=0) = gaps(gaps<=0)+L;
end

% -------------------------------------------------------------------------

function h = plotcirc(L,x,dt,h)
W = 0.05;  R = 1;
ncar = length(x);

theta = [0 : 2*pi/L : 2*pi];
xc = cos(theta);     yc = sin(theta);
xinner = (R-W/2)*xc; yinner = (R-W/2)*yc;
xouter = (R+W/2)*xc; youter = (R+W/2)*yc;

xi = [xinner(x);  xinner(x+1); xouter(x+1); xouter(x)];
yi = [yinner(x);  yinner(x+1); youter(x+1); youter(x)];
if nargin == 3
    color = randperm(ncar);
    h = fill(xi,yi, color); hold on
    plot(xinner,yinner, 'k', xouter,youter, 'k','linewidth',1.5)
    plot([xinner; xouter], [yinner; youter],'k','linewidth',1.5)
    axis image; axis((R+2*W)*[-1 1 -1 1]); axis off
else
    for i=1:ncar;  set(h(i),'xdata',xi(:,i),'ydata',yi(:,i));  end
end
pause(dt)