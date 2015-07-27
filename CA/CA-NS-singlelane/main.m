% This the main script implements the Nagel Schreckenberg cellular automata
% based traffic model. Flow vs density curves are displayed. 
%
% zhou lvwen: zhou.lv.wen@gmail.com
% 
density = 0:0.02:1;
roadlength = 100;
vmax = 5;
tmax = 200;
pbrak = 0;
flux = [];
vmean = [];

for rho = density
    [R, J, V] = ns(rho, pbrak, roadlength, tmax, 0, 0);
    flux = [flux; J];
    vmean = [vmean; V];
end

% ------------------------- density vs. volecity --------------------------
figure
plot(density, vmean,'k.','markersize',15);
hold on
plot(density,min(vmax, 1./density-1),'-r','linewidth',2)
ylim([0,5.55])
legend({'Cellular automata aproach', ...
        '$v(\rho) = \min\{v_{\max}, 1/\rho-1\}$'}, ...
        'interpreter','latex')
xlabel('density in vehicles/cell')
ylabel('velocity in cell/time')


% --------------------------- density vs. flux ----------------------------
figure
plot(density, flux,'k.','markersize',15);
hold on;
plot(density,min(density*vmax, 1-density),'-r','linewidth',2)
legend({'Cellular automata aproach', ...
        '$J(\rho) = \min\{\rho\cdot v_{\max}, 1-\rho\}$'}, ...
        'interpreter','latex')
xlabel('density in vehicles/cell')
ylabel('flux in vehicles/time')