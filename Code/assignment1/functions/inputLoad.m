function [input] = inputLoad
% inputLoad
% 
% Function to load the input structure for assignment 1 into workspace.
% 
% PROTOTYPE:
%  [input] = inputLoad
% 
% INPUT:
% 
% OUTPUT:
%  input    Input structure
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  02-11-2019: First version


% Assignment data & constrains:
% Given Planets:
input.id_d = 8;                               % Neptune
input.id_f = 4;                               % Mars
input.id_a = 3;                               % Earth

% Planet Parameters:
input.ksun = astroConstants(4);               % Sun gravity constant 

input.k_d = astroConstants(18);               % Neptune gravity constant
input.k_f = astroConstants(14);               % Mars gravity constant
input.k_a = astroConstants(13);               % Earth gravity constant

input.Rpl_d = astroConstants(28);             % Neptune mean radius
input.Rpl_f = astroConstants(24);             % Mars mean radius
input.Rpl_a = astroConstants(23);             % Earth mean radius

input.hatm_f = 250;                           % Height of Mars atmosphere [km]

input.TsynME = 778.77;           % Mars and Earth synodic period [julian days] 

% Time window
input.ED = date2mjd2000([2020 1 1 0 0 0]);    % Earliest Departure
input.LA = date2mjd2000([2060 1 1 0 0 0]);    % Latest Arrival

% Setup for gridSearch.m:

% TOF limits:
input.gridSearch.minTOF1 = 13.36*365;         % Minimum time of flight from Neptune to Mars [days]
input.gridSearch.maxTOF1 = 37.76*365;         % Maximum time of flight from Neptune to Mars [days]
input.gridSearch.minTOF2 = 24.13;             % Minimum time of flight from Mars to Earth [days]
input.gridSearch.maxTOF2 = 310.64;            % Maximum time of flight from Mars to Earth [days]

% Time steps:
input.gridSearch.dt_d = 165;                % Departure time [days]
input.gridSearch.dt_f = 1;                 % Flyby time [days]
input.gridSearch.dt_a = 1;                 % Arrival time [days]
 
% Setup for gaSearch.m:
% TOF limits:
input.gaSearch.minTOF1 = 13.36*365;           % Minimum time of flight from Neptune to Mars
input.gaSearch.maxTOF1 = 37.76*365;           % Maximum time of flight from Neptune to Mars
input.gaSearch.minTOF2 = 1*365;            % Minimum time of flight from Mars to Earth
input.gaSearch.maxTOF2 = 12*365;           % Maximum time of flight from Mars to Earth

% Parameters:
input.gaSearch.ngen = 300;                   % Number of generation
input.gaSearch.npop = 2000;                   % Population size
input.gaSearch.ConstraintTolerance = 1e-3;   % Constraint tolerance