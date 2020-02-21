function [par1,par2] = transferArcs(x,input)
% 
% transferArcs Transfer arcs parameters
% 
% Function to compute the keplerian elements of the Lambert arcs
% 
% PROTOTYPE:
%  [par1,par2] = transferArcs(x,input)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: x(1) = departure
%                                      x(2) = flyby
%                                      x(3) = arrival
%  input          Input structure
% 
% OUTPUT:
%  par1 [7,1]     Parameters for first arc 
%  par2 [7,1]     Parameters for second arc
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

t_d = x(1); t_f = x(2); t_a = x(3); % Initialize times

% Planets id numbers:
id_d = input.id_d;                              % Neptune
id_f = input.id_f;                              % Mars
id_a = input.id_a;                              % Earth

[kep_d,ksun] = uplanet(t_d, id_d);
[rr_d] = par2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, id_f);
[rr_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,ksun] = uplanet(t_a, id_a);
[rr_a] = par2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,~,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,~,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

[~,~,~,~,~,th1_i,kep1] = car2par(rr_d,vt1_i',input.ksun);
[~,~,~,~,~,th1_f] = car2par(rr_f,vt1_f',input.ksun);
[~,~,~,~,~,th2_i,kep2] = car2par(rr_f,vt2_i',input.ksun);
[~,~,~,~,~,th2_f] = car2par(rr_a,vt2_f',input.ksun);

par1 = [kep1(1:5);th1_i;th1_f];
par1(3:end) = rad2deg(par1(3:end));
par2 = [kep2(1:5);th2_i;th2_f];
par2(3:end) = rad2deg(par2(3:end));
