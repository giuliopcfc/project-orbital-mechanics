function dv = dvFun(x, input)
%
% dvFun Total deltaV of the mission
% 
% Function to compute the total deltaV of the mission given the times of
% departure, flyby and arrival. In this case, there is no
% constraint on the minimum value of radius of pericenter of the flyby
% manouvre (rp_min = 0), this constraint will be imposed with the function
% rpCon.m.
% 
% PROTOTYPE:
%  dv = dvFun(x)
% 
% INPUT:
%  x [3,1]        Times array:         x(1) = departure time [MJD2000]
%                                      x(2) = first transfer tof [days]
%                                      x(3) = second transfer tof [days]
% 
% OUTPUT:
%  dv [1]         Total deltaV of the mission [km/s]
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

% Initialize times:
x1 = cumsum(x);
t_d = x1(1); t_f = x1(2); t_a = x1(3); 

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

dv = NaN; % Initialize output

% Check: there must be no errors in lambertMR
if (err1==0)&&(err2==0)
    
    vinfm = vt1_f' - vv_f;
    vinfp = vt2_i' - vv_f;
    
    dvp = powGA(vinfm,vinfp,0,0, input);
    
    dv = norm(vv_d - vt1_i')+norm(vv_a - vt2_f') + dvp;
    
end
