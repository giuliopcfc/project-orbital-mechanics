function [dv_ga, x_ga] = gaSearch(input)
% gaSearch
% 
% Function to find the minimum delta velocity and the related time of
% departure and TOFs of the two transfer arcs using the genetic algorithm.
% It takes as input the boundaries of the whole time window, the restricted
% slots of times in which evaluate the population and where it will
% generate some tries in order to minimize the global velocity changes
% given by the engine when it performes manoeuvres.
% 
% PROTOTYPE:
%  [dv_ga, x_ga] = gaSearch(input)
%  
% INPUT:
%  input    Input structure
% 
% OUTPUT:
%  dv_ga   Minimum deltaV found [km/s]
%  x_ga[3,1]    MJD2000 times array for minimum deltaV: x_ga(1) = departure
%                                                       x_ga(2) = flyby
%                                                       x_ga(3) = arrival 
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  04-12-2019: First version
% 

% Setting time boundaries

% TOF limits:
minTOF1 = input.gaSearch.minTOF1;         % Minimum time of flight from Neptune to Mars
maxTOF1 = input.gaSearch.maxTOF1;         % Maximum time of flight from Neptune to Mars
minTOF2 = input.gaSearch.minTOF2;         % Minimum time of flight from Mars to Earth
maxTOF2 = input.gaSearch.maxTOF2;         % Maximum time of flight from Mars to Earth

% Time window
ED = input.ED;                               % Earliest Departure
LA = input.LA;                               % Latest Arrival

LB = [ED;                                    % Lower boundary
      minTOF1;
      minTOF2];
UB = [ED + input.TsynME;                     % Upper boundary
      maxTOF1
      maxTOF2];

%% Linear Costraint
% Building the matrix and the vector of known terms for composing the
% linear system which introduces the time boundaries in ga
A = [1 1 1];                                 
B = [LA];                                    

%% Genetic Algorithm
% After having set the parameters about how many tries evaluate and generate
% ga computes the absolute minimum delta velocity and it displays the times
ngen = input.gaSearch.ngen;                           % Number of generations
npop = input.gaSearch.npop;                           % Population size

options = optimoptions(@ga,'MaxGeneration',ngen,...
    'PopulationSize',npop,...
    'NonlinConAlgorithm','penalty',...
'ConstraintTolerance',input.gaSearch.ConstraintTolerance,...
'FunctionTolerance',0.01);

% Genetic Algorithm
[x_ga,dv_ga] = ga(@(x) dvFun(x, input),3,A,B,[],[],LB,UB,@(x) rpCon(x,input),options);

x_ga = cumsum(x_ga);                         % Converte x_ga to date array
