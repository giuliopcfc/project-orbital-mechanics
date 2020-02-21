function [dv,dv1,dv2,dv3,rp,dv_fb, vt1_f, vt2_i] = totdv(x, input)
%
% Function to compute the total deltaV of the mission given the times of
% departure, flyby and arrival with t_d<t_f<t_a. 
% 
% PROTOTYPE:
%  [dv,dv1,dv2,dv3,rp,dv_fb, vt1_f, vt2_i] = totdv(x, input)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: x(1) = departure
%                                      x(2) = flyby
%                                      x(3) = arrival
%  input          Input structure
% 
% OUTPUT:
%  dv [1]         Total deltaV of the mission [km/s]
%  dv1 [1]        First deltaV [km/s]
%  dv2 [1]        Second deltaV [km/s]
%  dv3 [1]        Third deltaV [km/s]
%  rp  [1]        Radius of pericenter of flyby manoevure [km]
%  dv_fb [1]      Total deltaV Flyby [km/s]
%  vt1_f [3,1]    Incoming velocity vector for flyby [km/s]
%  vt2_i [3,1]    Outgoing velocity vector for flyby [km/s]
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  09-12-2019: First version
% 
               
t_d = x(1); t_f = x(2); t_a = x(3); % Initialize times

% Load constants:

% Planets id numbers:
id_d = input.id_d;                              % Neptune
id_f = input.id_f;                              % Mars
id_a = input.id_a;                              % Earth

[kep_d,ksun] = uplanet(t_d, id_d);
[rr_d,vv_d] = par2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, id_f);
[rr_f,vv_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,ksun] = uplanet(t_a, id_a);
[rr_a,vv_a] = par2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,err1,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,err2,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

dv = NaN; dv1 = NaN; dv2 = NaN; dv3 = NaN; rp = NaN; % Initialize output

% Check: there must be no errors in lambertMR
if (err1==0)&&(err2==0)
    
    % Mars Reference frame
    vinfm = vt1_f' - vv_f;
    vinfp = vt2_i' - vv_f;
    
    [dvp,~,rp] = powGA(vinfm,vinfp,input.Rpl_f, input.hatm_f, input);
    
    dv = norm(vv_d - vt1_i')+norm(vv_a - vt2_f') + dvp;
    dv1 = norm(vv_d - vt1_i');
    dv2 = dvp;
    dv3 = norm(vv_a - vt2_f');
    dv_fb = norm(vinfm-vinfp);
    
end

end