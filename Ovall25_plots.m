%% ------------------------------------------------------------------------
% Script: LeConte Glacier Vertical Section Analysis
% Purpose: Load T, S, U, V, W fields from NetCDF outputs and produce plots
%          of vertical sections and plume transport diagnostics.
% -------------------------------------------------------------------------

%% --- File setup ----------------------------------------------------------
filepre = 'leconte_vert1_d50w_nam2';   % base filename prefix for this run

% Scalars
T = ncread([filepre '_t.nc'],'t');     % Temperature
S = ncread([filepre '_s.nc'],'s');     % Salinity

% Velocity components
u = ncread([filepre '_u.nc'],'u');     % x-velocity
v = ncread([filepre '_v.nc'],'v');     % y-velocity
w = ncread([filepre '_w.nc'],'w');     % vertical velocity

% Grid variables (load from one representative file)
xw   = ncread([filepre '_v.nc'],'xw'); % x-coordinates for velocity grid
yu   = ncread([filepre '_v.nc'],'yv'); % y-coordinates
zw   = ncread([filepre '_v.nc'],'zu'); % vertical coordinates
time = ncread([filepre '_v.nc'],'time');

% Grid dimensions
Ny = size(u,2);    % number of y-points
midY = Ny/2;       % index of the midline in y-direction

%% --- Temperature profile near the wall ----------------------------------
figure
plot(squeeze(T(8,midY,17:end)), zw(17:end));
xlim([5 8])
xlabel('Temperature (°C)')
ylabel('z (m)')
title('Temperature at x = 0.37 m from wall')
grid on

%% --- Vertical sections of velocity --------------------------------------
% v at a vertical slice, x = 5.6 m
figure
pcolor(yu, zw, squeeze(v(8,:,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
title('v at x = 5.6 m')

% v at midline
figure
pcolor(xw, zw, squeeze(v(:,midY,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
xlim([0 100]), title('v at midline')

% u at midline
figure
pcolor(xw, zw, squeeze(u(:,midY,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
xlim([0 100]), title('u at midline')

% u at x = 6 m
figure
pcolor(yu, zw, squeeze(u(9,:,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
title('u at x = 6 m')

% w at x = 5.6 m
figure
pcolor(yu, zw, squeeze(w(8,:,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
title('w at x = 5.6 m')

% w at x = 7.1 m
figure
pcolor(yu, zw, squeeze(w(10,:,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
title('w at x = 7.1 m')

% w at midline with contours
figure
pcolor(xw, zw, squeeze(w(:,midY,:))'), shading interp
cmocean('balance'), caxis([-.1 .1])
title('w at midline')
hold on
contour(xw, zw, squeeze(w(:,midY,:))',[0.01:.1:2],'k')

%% --- Transport calculations ---------------------------------------------
% Full domain transport (scaled by grid spacing ~0.75 m)
transport_total = squeeze(nansum(nansum(w))) * (0.75^2);

% Focused transport near plume (subset indices)
transport_plume = squeeze(nansum(nansum(w(1:50, midY-60:midY+60,:)))) * (0.75^2);

% Plot plume transport vs depth
figure
plot(transport_plume, zw)
title('Vertical transport within plume')
xlabel('Transport (m^3/s per level)')
ylabel('z (m)')

%% --- Vertical profiles at selected depths -------------------------------
figure
plot(xw, w(:,midY,67),  'DisplayName','z = 50 m'); hold on
plot(xw, w(:,midY,133), 'DisplayName','z = 100 m')
plot(xw, w(:,midY,173), 'DisplayName','z = 150 m')
legend show
xlim([5 50])
title('Vertical velocity profiles at midline')

%% --- Crossflow transports (optional) -----------------------------------
ytransport = squeeze(nansum(v(:,435,:))) * 0.75;   % along-y transport
xtransport = squeeze(nansum(u(50,:,:),2)) * 0.75;  % along-x transport

