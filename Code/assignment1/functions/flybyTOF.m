function [TOF_f, data_f] = flybyTOF(x, input)
% 
% flybyTOF time of flight of the flyby
% 
% Function to compute the time of flight of the flyby given the date of
% departure, flyby and arrival and other parameters of the flyby.
% 
% PROTOTYPE:
%  [TOF_f] = flybyTOF(x, input)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: x(1) = departure
%                                      x(2) = flyby
%                                      x(3) = arrival
%  input          Input structure
% 
% OUTPUT:
%  TOF_f [1]      Time duration of flyby [s]
%  data_f         Structure of other relevant data for flyby
%
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  01-12-2019: First version
%

t_f = x(2); % Flyby date [MJD]

% Mars parameters at flyby:
[kep_f,ksun] = uplanet(t_f, input.id_f); % Keplerian elements
[rr_f,vv_f] = par2car(kep_f(1),kep_f(2),...
    kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun); % Cartesian elements

[~,~,~,~,~,~, VFm, VFp] = totdv(x, input); % Get incoming and outgoing velocity vectors

vinfm = VFm' - vv_f; % Excess velocity vector of the incoming hyperbola
vinfp = VFp' - vv_f; % Excess velocity vector of the outgoing hyperbola

[~,~, rp, em, ep, am, ap, vpp, vpm] =...
    powGA(vinfm,vinfp,input.Rpl_f, input.hatm_f, input); % Compute powered gravity assist

% Planets' gravity constants:
ksun = input.ksun;
kmars = input.k_f;

r_SOI = norm(rr_f)*(kmars/ksun)^(2/5); % Radius of Mars SOI [km]

% Apply hyperbola time law:

% Incoming branch:
fm = acos((am*(1-em^2)-r_SOI)/em/r_SOI); % True anomaly at SOI [rad]
Fm = 2*atanh(tan(fm/2)/sqrt((1+em)/(em-1))); % Eccentric anomaly at SOI [rad]
TOFm = sqrt(-(am^3)/kmars)*(em*sinh(real(Fm))-real(Fm)); % TOF of incoming branch [s]

% Outgoing branch:
fp = acos((ap*(1-ep^2)-r_SOI)/ep/r_SOI); % True anomaly at SOI [rad]
Fp = 2*atanh(tan(fp/2)/sqrt((1+ep)/(ep-1))); % Eccentric anomaly at SOI [rad]
TOFp = sqrt(-(ap^3)/kmars)*(ep*sinh(real(Fp))-real(Fp)); % TOF of outgoing branch [s]

TOF_f = TOFm + TOFp; % Time of flight of flyby [s]

% Save flyby data:
data_f.r_SOI = r_SOI;
data_f.em = em;
data_f.ep = ep;
data_f.vinfm = norm(vinfm);
data_f.vinfp = norm(vinfp);
data_f.VFm = VFm;
data_f.VFp = VFp;
data_f.TOFm = TOFm;
data_f.TOFp = TOFp;
data_f.vpm = vpm;
data_f.vpp = vpp;
data_f.rp = rp;
