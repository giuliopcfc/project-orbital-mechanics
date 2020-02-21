function [] = plotFlyby(x,name,input)
% 
% Function to plot the flyby given the times of departure, flyby and 
% arrival.
% 
% PROTOTYPE:
%  plotFlyby(x)
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

[~, data_f] = flybyTOF(x, input);
    
vpp = data_f.vpp*[0;1;0];
vpm = data_f.vpm*[0;1;0];
rp = data_f.rp*[1;0;0];
TOFm = data_f.TOFm;
TOFp = data_f.TOFp;

y0m = [rp; vpm];
[~,ym] = twoBodyInt([0 -TOFm],y0m,input.k_f);

y0p = [rp; vpp];
[~,yp] = twoBodyInt([0 TOFp],y0p,input.k_f);

h = figure;
hold on
planet3d(4, [0 0 0], h, 1)
axis equal
grid on
plot3(ym(:,1),ym(:,2),ym(:,3))
plot3(yp(:,1),yp(:,2),yp(:,3))
planet3d(4, [0 0 0], h, 1)


title(name)

