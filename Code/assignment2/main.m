%% main.m Main script for assignment 2
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici
%  Luca Rizzieri
%  Davide Sisana
% 
% VERSIONS:
%  03-02-2019: First version
%

%% Initialization:

addpath([pwd,'\functions'])

clear;clc;close all;

%% Load input data and constants:

input = inputLoad; % Input structure

%% Ground tracks:

times = [input.T, 24*3600, 10*24*3600]; % Required times of integration [s]

input.tspan = [0:0.5:times(end)]; % Set time span array for orbit integration

% Find indexes of input.tspan that corresponds to required times:
[~,i1] = min(abs(input.tspan - input.T));
[~,i2] = min(abs(input.tspan - 24*3600));

indexes = [i1;
    i2;
    length(input.tspan)];

% Orbit propagation:

% Unperturbed:
input.type = 0;

[lonNaN0,latNaN0] = GT(input);

% With J2 effect:
input.type = 1;

[lonNaN1,latNaN1] = GT(input);

% Plot ground tracks:

strings = {'1 orbit','1 day','10 days'}; % Plot titles

for i = 1:3

figure,
hold on
plot(lonNaN0(1:indexes(i)),latNaN0(1:indexes(i)),'g','linewidth',1.5)
plot(lonNaN0(1),latNaN0(1),'ro','linewidth',3)
plot(lonNaN0(indexes(i)),latNaN0(indexes(i)),'rx','linewidth',3)
EarthGroundTrack
title(['Unperturbed GT ', strings{i}])
hold off

figure,
hold on
plot(lonNaN1(1:indexes(i)),latNaN1(1:indexes(i)),'y','linewidth',1.5)
plot(lonNaN1(1),latNaN1(1),'ro','linewidth',3)
plot(lonNaN1(indexes(i)),latNaN1(indexes(i)),'rx','linewidth',3)
EarthGroundTrack
title(['Perturbed GT ', strings{i}])
hold off

end

%% Repeating ground tracks:

% Get semimajor axis values for repeating ground tracks:
[a_pert,a_unpert] = repGTJ2(input)

% Time span array from 0 to m Earth rotations:
input.tspan = [0:0.5:input.m*2*pi/norm(input.w_E)];

% Orbit propagation:

% Unperturbed:
input.type = 0;

% Update initial state array for GPE method:
input.Y0GPE(1) = a_unpert;

[lonNaN0_rep,latNaN0_rep] = GT(input);

% With J2 effect:
input.type = 1;

% Update initial state array for GPE method:
input.Y0GPE(1) = a_pert;

[lonNaN1_rep,latNaN1_rep] = GT(input);

% Plot repeating ground tracks:

figure,
hold on
plot(lonNaN0_rep,latNaN0_rep,'g','linewidth',1.5)
plot(lonNaN0_rep(1),latNaN0_rep(1),'ro','linewidth',3)
plot(lonNaN0_rep(end),latNaN0_rep(end),'rx','linewidth',3)
EarthGroundTrack
title(['Unperturbed Repeating GT'])
hold off

figure,
hold on
plot(lonNaN1_rep,latNaN1_rep,'y','linewidth',1.5)
plot(lonNaN1_rep(1),latNaN1_rep(1),'ro','linewidth',3)
plot(lonNaN1_rep(end),latNaN1_rep(end),'rx','linewidth',3)
EarthGroundTrack
title(['Perturbed Repeating GT'])
hold off

%% Accuracy of the two integration methods:
input = inputLoad;

input.type = 3;

input.tspan = [0:500:100*input.T]; % Time span of one year

% Propagation with cartesian method:
input.method = 'cart';

[t,kep_c,cart_c] = orbitIntegration(input);

% Propagation with GPE method:
input.method = 'GPE';

[t,kep_GPE,cart_GPE] = orbitIntegration(input);

% Wrap from 0 to 2*pi and convert angles to degrees:
kep_c(:,3:end) = rad2deg(wrapTo2Pi(kep_c(:,3:end)));
kep_GPE(:,3:end) = rad2deg(wrapTo2Pi(kep_GPE(:,3:end)));


% Error on keplerian elements:
err_kep = kep_GPE - kep_c; 
err_kep(:,3:end) = wrapTo180(err_kep(:,3:end));
err_kep = max(abs(err_kep))

%% Comparison between computational times required by the two methods:

input = inputLoad;

input.type = 3; % Set perturbation type

norbits = 100; % Number of orbits

tspan_max = norbits*input.T; % Maximum 

npoints = norbits*[200:20:1e4]; % Array of number of time steps for each integration

t_cart = []; % Array of computational times for cartesian method
t_GPE = []; % Array of computational times for GPE method

for i = 1:length(npoints)
    
    input.tspan = linspace(0,tspan_max,npoints(i));
   
    % Cartesian:
    tic;
    input.method = 'cart';
    orbitIntegration(input);
    t_cart = [t_cart,toc]
    
        % Gauss planetary equations:
    tic;
    input.method = 'GPE';
    orbitIntegration(input);
    t_GPE = [t_GPE, toc]
    
end

figure,
i_cut = 200; % Cut the first 200 samples
semilogy(npoints(i_cut:end),t_cart(i_cut:end),npoints(i_cut:end),t_GPE(i_cut:end)), 
legend('Cartesian','GPE'),
title('Computational time'),
xlabel('Number of time steps'),
ylabel('Time [s]'),
grid on
%% Time evolution of keplerian elements:

input.type = 3;

input.tspan = [0:500:10*365*24*3600]; % Time span of one year

% Propagation with GPE method:
input.method = 'GPE';

[t,kep_GPE,cart_GPE] = orbitIntegration(input);

kep_GPE(:,3:end) = rad2deg(wrapTo2Pi(kep_GPE(:,3:end)));

% Plot keplerian elements history:

labels = {'Semi major axis [km]','Eccentricity [-]','Inclination [deg]','RAAN [deg]',...
    'Anomaly of pericenter [deg]','True anomaly [deg]'};
figure,
hold on
for i = 1:6
    
subplot(3,2,i)
plot(t/365/24/3600,kep_GPE(:,i))
ylabel(labels{i})
xlabel('Time [years]')

end
hold off

%% Orbit evolution:

input = inputLoad;

N_orbits = 52; % Number of orbits to plot

% Orbit integration:
input.tspan = linspace(0,365*24*3600,N_orbits); 
input.method = 'GPE';
[t,kep] = orbitIntegration(input);

% Plot orbit evolution:
h = figure;
hold on

myCMap = parula(length(t));
set(gcf,'DefaultAxesColorOrder', myCMap);
colormap(myCMap);
hcb = colorbar;
set(get(hcb,'Title'),'String','\Deltat [days]');
caxis([0 t(end)]/3600/24)

for i = 1:N_orbits
    
    orbCol = myCMap(i,:);
    
    [X,Y,Z] = plotOrbit(kep(i,1),kep(i,2),kep(i,3),kep(i,4),kep(i,5),...
        0,2*pi,0.01,input.k_E);
    
    plot3(X,Y,Z,'color',orbCol);
    
    plot3([X(1) X(round(end/2))],[Y(1) Y(round(end/2))],...
        [Z(1) Z(round(end/2))],'color',orbCol,...
        'linestyle','--');
    xlabel('x [km]'),ylabel('y [km]'),zlabel('z [km]'),
    
    title('Orbit evolution')
    
end
planet3d(3, [0,0,0], h, 1);
hold off
axis equal

%% Frequency analysis:

% Orbit integration:

dt = 500; % Sampling time [s]
T = 10*365*24*3600; % Time of observation [s]

input.type = 3;
input.method = 'GPE';
input.tspan = [0:dt:T];

% Check if input.tspan has an even number of elements, if not delete the
% last element.
if mod(length(input.tspan),2) 
    input.tspan(end) = [];
end

[t,kep,cart] = orbitIntegration(input); % Integration
kep(:,3:6) = rad2deg(kep(:,3:6)); % Convert angles from radians to degrees

N = length(t); % Number of samples
fs = 1/dt; % Sampling frequency [Hz]

% Apply fft:
y = fft(kep,N,1); 
ynorm = abs(y/N);

% Keplerian elements Fourier transform (single-sided spectrum):
f_kep = ynorm(1:N/2+1,:);
f_kep(2:end-1,:) = 2*f_kep(2:end-1,:);

f_ax = [0:N/2]*fs/N; % Frequency axis [Hz]

% Plot Fourier transforms without the frequency 0 Hz:

figure,

ylabels = {'Semi major axis [km]','Eccentricity [-]','Inclination [deg]',...
    'RAAN [deg]', 'Anomaly of pericenter [deg]'};

for i = 1:5
    
    subplot(3,2,i),
    stem(f_ax(2:end),f_kep(2:end,i))
    xlabel('Frequency [Hz]'), ylabel(ylabels{i})
end

%% Filtering:

fco_td = 1e-8*[1 1 1 6 6]; % Cut-off frequencies array for time domain filter [Hz]
fco_fd = 1e-6*[1 1 1 4 4]; % Cut-off frequencies array for frequency domain filter [Hz]

% Filtering with centered moving mean:

kep_td = kep; % Initialize matrix of filtered elements

for i = 1:5
    
    k = fs/fco_td(i)/2; % Number of elements on which perform the mean
    
    kep_td(:,i) = movmean(kep(:,i),k);
    
end

% Filtering with low pass filter in frequency domain:

y_fd = y; % Initialize filtered Fourier transform

for i = 1:5
    
    % Find index of f_ax that corresponds to the cut off frequency:
    [~,index] = min(abs(f_ax - fco_fd(i))); 
    
    % Set to zero the values of Fourier transform that corresponds to
    % frequencies higher than the cut off frequency:
    y_fd(index:end,i) = 0;
    
end

kep_fd = ifft(y_fd,'symmetric'); % Inverse transform of filtered Fourier transform

% Plot filtered vs unfiltered:
figure,

ylabels = {'Semi major axis [km]','Eccentricity [-]','Inclination [deg]',...
    'RAAN [deg]', 'Anomaly of pericenter [deg]'};

kep(:,3:5) = wrapTo360(kep(:,3:5));
kep_fd(:,3:5) = wrapTo360(kep_fd(:,3:5));
kep_td(:,3:5) = wrapTo360(kep_td(:,3:5));

for i = 1:5
    
    subplot(3,2,i),
    hold on,
    plot(t/24/3600,kep(:,i))
    plot(t/24/3600,kep_fd(:,i),'r','linewidth',1.5)
    plot(t/24/3600,kep_td(:,i),'g','linewidth',1.5)
    legend('Unfiltered','Frequency domain filter','Time domain filter')
    ylabel(ylabels{i})
    xlabel('Time [days]')
    hold off
end

%% Real data comparison:

% Load data from horizon_results.csv:
csv_mat = csvread('real_data\horizons_results.csv');
t = (csv_mat(:,1) - csv_mat(1,1))*24*3600; % Time axis

% Extract keplerian elements:
a_real = csv_mat(:,11); % Semi major axis [km]
e_real = csv_mat(:,2); % Eccentricity
i_real = csv_mat(:,4); % Inclination [deg]
OM_real = csv_mat(:,5); % RAAN [deg]
om_real = csv_mat(:,6); % Anomaly of pericenter [deg]
f_real = csv_mat(:,10); % True anomaly [deg]

% Set initial conditions for integration:
input.a0 = a_real(1);
input.e0 = e_real(1);
input.i0 = i_real(1)*pi/180;
input.OM0 = OM_real(1)*pi/180;
input.om0 = om_real(1)*pi/180;
input.f0 = f_real(1)*pi/180;

input.Y0GPE = [input.a0;input.e0;input.i0;input.OM0;...
    input.om0;input.f0]; % Initial state array for GPE

% Real orbit propagation with Am = 0.0054:
input.tspan = t;
input.method = 'GPE';
input.type = 3;
input.Am = 0.0054;

[t,kep] = orbitIntegration(input);

a1 = kep(:,1); 
e1 = kep(:,2);
i1 = wrapTo2Pi(kep(:,3))*180/pi;
OM1 = wrapTo2Pi(kep(:,4))*180/pi;
om1 = wrapTo2Pi(kep(:,5))*180/pi;

% Comparison:
figure,
subplot(3,2,1)
plot(t/24/3600,a1,t/24/3600,a_real)
xlabel('Time [days]'), ylabel('a [km]')
legend('MATLAB model A/m = 0.0054','NASA HORIZONS')
grid on

subplot(3,2,2)
plot(t/24/3600,e1,t/24/3600,e_real)
xlabel('Time [days]'), ylabel('e [-]')
legend('MATLAB model A/m = 0.0054','NASA HORIZONS')
grid on

subplot(3,2,3)
plot(t/24/3600,i1,t/24/3600,i_real)
xlabel('Time [days]'), ylabel('i [deg]')
legend('MATLAB model A/m = 0.0054','NASA HORIZONS')
grid on

subplot(3,2,4)
plot(t/24/3600,OM1,t/24/3600,OM_real)
xlabel('Time [days]'), ylabel('\Omega [deg]')
legend('MATLAB model A/m = 0.0054','NASA HORIZONS')
grid on

subplot(3,2,5)
plot(t/24/3600,om1,t/24/3600,om_real)
xlabel('Time [days]'), ylabel('\omega [deg]')
legend('MATLAB model A/m = 0.0054','NASA HORIZONS')
grid on

