function [rho, flux, vmean, Nacdnts] = nsacdnt(rho, p, L, tmax, animation, spacetime)
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
%   4. Random speed up.  vn → vn + 1 presuming vn < vmax
%
%   5. Vehicle movement. The vehicles are deterministically moved by their 
%      velocities, xn → xn + vn .
%
% USAGE: flux = nsacdnt(rho, p, L, tmax, isdraw)
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
rand('seed',0)
if nargin == 0; 
    rho = 0.25; p = 0.25; L = 100; tmax = 200; pacdnt = 0.001;
    animation = true; spacetime = true;
end

vmax = 5;                  % maximun speed
% place a distribution with density
ncar = round(L*rho);
rho = ncar/L;

x = sort(randperm(L, ncar));
v = vmax * ones(1,ncar);   % start everyone initially at vmax
isacdnt = zeros(size(x));

if animation; h = plotcirc(L,x,isacdnt,0.1); end

flux = 0;                  % number of cars that pass through the end
vmean = 0;
Nacdnts = 0;
road = zeros(tmax, L);
for t = 1:tmax
    % acceleration
    v = min(v+1, vmax);

    %collision prevention
    gaps = gaplength(x,L); % determine the space vehicles have to move

    v = min(v, gaps-1);
    
    % random speed drops
    vdrops = ( rand(1,ncar)<p & v>0);
    v = v-vdrops;
    
    % random speed add
    vadds = ( rand(1,ncar)<pacdnt & v>=0);
    v = min(v+vadds,vmax);
    
    
    % update the position
    x(v>0) = x(v>0) + v(v>0);
    passed = x>L;          % cars passed at time r
    x(passed) = x(passed) - L;% periodic boundary conditions
    
    [isacdnt, v, nacdnt] =  accident(isacdnt, v, x);
    Nacdnts = Nacdnts + nacdnt;
    if t>tmax/2
        flux = flux + sum(v(v>0)/L); %flux = flux + sum(passed); 
        vmean = vmean + mean(v(v>0));
    end
    road(t,x) = 1;
    if animation; [h,ha] = plotcirc(L,x,isacdnt,0.1,h); delete(ha);end
    
end
flux = flux/(tmax/2);
vmean = vmean/(tmax/2);

if spacetime; figure;imagesc(road);colormap([1,1,1;0,0,0]);axis image; end

% -------------------------------------------------------------------------

function [isacdnt, v, nacdnt] =  accident(isacdnt, v, x)
for xi =  unique(x);
    tot = find(x==xi);             % total
    old = find(x==xi&isacdnt==1);  % old
    new = find(x==xi&isacdnt==0);  % new
    if length(tot)>=2
        if isempty(old) && ~isempty(new)
            vmin = -1;
            v(new) = vmin - 5*[0:length(new)-1];
        else
            vmin = min(-1,min(v(old)));
            v(new) = vmin - 5*[1:length(new)];
        end
       isacdnt(tot) = 1;
    end
end
isacdnt(v>0) = 0;
nacdnt = sum(isacdnt);

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

function [h, ha] = plotcirc(L,x,isacdnt,dt,h)
W = 0.05;  R = 1;
ncar = length(x);

theta = [0 : 2*pi/L : 2*pi];
xc = cos(theta);     yc = sin(theta);
xinner = (R-W/2)*xc; yinner = (R-W/2)*yc;
xouter = (R+W/2)*xc; youter = (R+W/2)*yc;

xi = [xinner(x);  xinner(x+1); xouter(x+1); xouter(x)];
yi = [yinner(x);  yinner(x+1); youter(x+1); youter(x)];

x = x(isacdnt==1);
xa = [xinner(x);  xinner(x+1); xouter(x+1); xouter(x)];
ya = [yinner(x);  yinner(x+1); youter(x+1); youter(x)];
if nargin == 4
    color = randperm(ncar);
    h = fill(xi,yi, color); hold on
    plot(xinner,yinner, 'k', xouter,youter, 'k','linewidth',1.5)
    plot([xinner; xouter], [yinner; youter],'k','linewidth',1.5)
    axis image; axis((R+2*W)*[-1 1 -1 1]); axis off
else
    for i=1:ncar;  set(h(i),'xdata',xi(:,i),'ydata',yi(:,i));  end
    ha = plot(xa,ya,'r','linewidth',3);
end
pause(dt)
