function [dy] = twoBodyOde(~,y,mu)
% twoBodyOde Two-body problem ode function
% 
% Input function for the ode integration of the two-body problem.
% 
% PROTOTYPE:
%  [dy] = twoBodyOde(t,y,mu)
%  
% INPUT:
%  y  [6,1]     Stete vector at a certain time [km] and [km/s]
%  mu  [1]      Gravitational parameter of the primary body [km^3/s^2]
% 
% OUTPUT:
%  dy  [6,1]    Time derivative of the state vector [km/s] and [km/s^2]
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  02-11-2019: First version
%
dy = zeros(6,1);
r = norm(y(1:3));
dy(1:3) = y(4:6);

dy(4:6) = -mu/r^3*y(1:3);
end

