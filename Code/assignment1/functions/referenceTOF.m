function [TOFH_outer, TOFH_inner, TPAR_outer, TPAR_inner] = referenceTOF(input)
% referenceTOF
% 
% Function to compute the reference time of flight for the two transfer
% arcs.
% 
% PROTOTYPE:
%  [TOFH_outer, TOFH_inner, TPAR_outer, TPAR_inner] = referenceTOF(input)
% 
% INPUT:
%  input          Input structure
% 
% OUTPUT:
%  TOFH_outer     TOF for Hohmann transfer between Neptune and Mars [years]  
%  TOFH_inner     TOF for Hohmann transfer between Earth and Mars [days]
%  TPAR_outer     TOF for parabolic transfer between Neptune and Mars [years]
%  TPAR_inner     TOF for parabolic transfer between Earth and Mars [days]
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  02-11-2019: First version
%% 

% Choose a reference time for keplerian elements:
ref_time = date2mjd2000([2020 1 1 0 0 0]);

% Planets id numbers:
id_d = input.id_d;
id_f = input.id_f;
id_a = input.id_a;

% Get keplerian elements at reference time:
[kepN,ksun] = uplanet(ref_time, id_d);
[kepM] = uplanet(ref_time, id_f);
[kepE] = uplanet(ref_time, id_a);

% Planets position: ideal case of circular orbits and an angle of 180
% degrees between planets.
RN = norm(kepN(1))*[-1 0 0]';                            % Neptune position

% Hohmann transfer between outer planets ( Mars and Neptune ).
aH_outer = 0.5*(kepN(1)+kepM(1));                        % Semi major axis

% The transfer must be close to 360 degrees. We set a small angle a in
% order to avoid errors in lambertMR.m.
a = 0.0001;                                              % Small angle 
RM = kepM(1)*[cos(a) sin(a) 0]';                         % Mars position

% Get parabolic tof from lambertMR.m:
[~,~,~,~,~,~,TPAR_outer] = lambertMR(RN,RM,100,ksun,0,0,0);

% Planets position: ideal case of circular orbits and an angle of 180
% degrees between planets.
RM = norm(kepM(1))*[1 0 0]';                             % Mars position

% Hohmann transfer between inner planets ( Mars and Earth).
aH_inner = 0.5*(kepM(1)+kepE(1));                        % Semi major axis

%  Parabolic transfer between inner planets ( Mars and Earth ):
RE = kepE(1)*[cos(a) sin(a) 0]';

[~,~,~,~,~,~,TPAR_inner] = lambertMR(RM,RE,100,ksun,0,0,0);

TOFH_outer = pi*sqrt(aH_outer^3/ksun)/3600/24/365;       % [years]
TPAR_outer = TPAR_outer/3600/24/365;                     % [years]
TOFH_inner = pi*sqrt(aH_inner^3/ksun)/3600/24;           % [days]
TPAR_inner = TPAR_inner/3600/24;                         % [days]