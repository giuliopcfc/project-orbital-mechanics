%% main.m Main script for assignment 1
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici
%  Luca Rizzieri
%  Davide Sisana
% 
% VERSIONS:
%  04-02-2019: First version
%

%% Initialization:

folder = 'assignment1';
string = erase(pwd,folder);
addpath([pwd,'\functions'])
addpath([pwd,'\functions\timeConversion\time'])
clear;clc;close all;

%% Load input data and constants:

input = inputLoad; % Input structure

%% Reference quantities:

% Reference times of flight:
[TOFH_outer, TOFH_inner, TPAR_outer, TPAR_inner] = referenceTOF(input)

% DeltaV of Hohmann transfer from Neptune to Earth:
dv_h = hohmannNE(input)

%% Grid search:

disp('Grid search running..')
tic

[dv_gs, x_gs] = gridSearch(input)

toc

% Plot:
name = 'Transfer arcs from grid search result';
plotTransfer(x_gs,name,input);

%% Genetic Algorithm:

N_runs = 5; % Number of runs
% To check algorithm convergence set N_runs > 1.

% Save results for each run:
x_runs = [];
dv_runs = [];

disp(['ga search with ',num2str(N_runs),' runs running..']);
tic

for i = 1:N_runs

[dv_ga, x_ga] = gaSearch(input);

dv_runs = [dv_runs; dv_ga]
x_runs = [x_runs; x_ga];

end

toc

% Convergence check if N_runs > 1:
if N_runs > 1
    
    max_difference = max(dv_runs) - min(dv_runs)
    
end

% Select minimum deltaV solution:

[~,index] = min(dv_runs);
dv_ga = dv_runs(index);
x_ga = x_runs(index,:);

% Plot:
name = 'Transfer arcs from ga result';
plotTransfer(x_ga,name,input);

%% Transfer arcs characterization:

[par1,par2] = transferArcs(x_ga,input);

%% Flyby analysis:

% Plot:
name = 'Flyby from ga result';
plotFlyby(x_ga,name,input);

% Time duration:
[TOF_f, data_f] = flybyTOF(x_ga, input);

%% Pork chop plots:

% Neptune-Mars transfer:
% Set time arrays:
t_d1 = linspace(date2mjd2000([2020 1 1 0 0 0]),...
    date2mjd2000([2060 1 1 0 0 0]),...
    500);                                       % Departure time
t_f1 = linspace(date2mjd2000([2020 1 1 0 0 0]),...
    date2mjd2000([2060 1 1 0 0 0]),...
    500);                                       % Flyby time
PC1(t_d1,t_f1,x_ga,x_gs, input);

% Mars-Earth transfer:
% Set time arrays:
t_f2 = linspace(date2mjd2000([2020 1 1 0 0 0]),...
    date2mjd2000([2060 1 1 0 0 0]),...
    500);                                       % Flyby time
t_a2 = linspace(date2mjd2000([2020 1 1  0 0 0]),...
    date2mjd2000([2060 1 1 0 0 0]),...
    500);                                       % Arrival time
PC2(t_f2,t_a2,x_ga,x_gs, input); 

