function dv = hohmannNE(input)
% hohmannNE Hohmann transfer from Neptune to Earth
% 
% Function to compute the cost of Hohmann transfer from Neptune to Earth
% 
% PROTOTYPE:
%  dv = hohmannNE(input)
% 
% INPUT:
%  input    Input structure
% 
% OUTPUT:
%  dv       Cost of the transfer [km/s]
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  02-11-2019: First version

ref_time = input.ED; % Reference time for keplerian elements

[kepN,ksun] = uplanet(ref_time, input.id_d); % Neptune keplerian elements

[kepE,ksun] = uplanet(ref_time, input.id_a); % Earth keplerian elements

vN = sqrt(ksun/kepN(1)); % Neptune velocity
vE = sqrt(ksun/kepE(1)); % Earth velocity

% Transfer orbit:

rp = kepE(1); % Radius of pericenter
ra = kepN(1); % Radius of apocenter

e = (ra - rp)/(ra + rp); % Eccentricity
p = 2*ra*rp/(ra+rp); % Semilatus rectus

vp = sqrt(ksun/p)*(1+e); % Velocity at pericenter
va = sqrt(ksun/p)*(1-e); % Velocity at apocenter

% Cost of manoeuvre:

dv = abs(vN - va) + abs(vE - vp);