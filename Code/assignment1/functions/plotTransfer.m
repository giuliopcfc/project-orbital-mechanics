function plotTransfer(x,name,input)
% 
% Function to plot the two transfer arcs given the times of
% departure, flyby and arrival.
% 
% PROTOTYPE:
%  plotTransfer(x)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: x(1) = departure
%                                      x(2) = flyby
%                                      x(3) = arrival
%  name           Figure title
%  input          Input structure
% 
% OUTPUT:
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

% Load constants:

% Planets id numbers:
id_d = input.id_d;                              % Neptune
id_f = input.id_f;                              % Mars
id_a = input.id_a;                              % Earth

% Transfer Arcs:

[kep_d,ksun] = uplanet(t_d, id_d);
[rr_d,vv_d] = par2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, id_f);
[rr_f,vv_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,ksun] = uplanet(t_a, id_a);
[rr_a,vv_a] = par2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,err1,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,err2,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

% Integration of transfer arcs:

[~,y1] = twoBodyInt([0 (t_f-t_d)*24*3600],[rr_d;vt1_i'],ksun);

[~,y2] = twoBodyInt([0 (t_a-t_f)*24*3600],[rr_f;vt2_i'],ksun);

% Planets orbits:

[X_d,Y_d,Z_d] = plotOrbit(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),0,2*pi,0.01,ksun);

[X_f,Y_f,Z_f] = plotOrbit(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),0,2*pi,0.01,ksun);

[X_a,Y_a,Z_a] = plotOrbit(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),0,2*pi,0.01,ksun);

% Plot

h = figure;
hold on

plot3(X_d,Y_d,Z_d,'r--')
plot3(X_f,Y_f,Z_f,'r--')
plot3(X_a,Y_a,Z_a,'r--')

plot3(y1(:,1),y1(:,2),y1(:,3))
plot3(y2(:,1),y2(:,2),y2(:,3))

% Planets and Sun:
planet3d(8, rr_d, h, 10000)
planet3d(4, rr_f, h, 10000)
planet3d(3, rr_a, h, 10000)
planet3d(0, [0;0;0], h, 100)

title(name)
grid on

axis equal
hold off

