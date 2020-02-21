function [dvp, delta, rp, em, ep, am, ap, vpm, vpp, deltam, deltap] = powGA(vinfm,vinfp,Rpl,hatm_f,input)
% powGA powered gravity assist
% 
% Function to compute the powered gravity assist.
% 
% PROTOTYPE:
%  [dvp, delta, rp, em, ep, am, ap, vpm, vpp, deltam, deltap] = powGA(vinfm,vinfp,Rpl,hatm_f,input)
%  
% INPUT:
%  vinfm [3]    Incoming relative velocity to planet vector [km/s]
%  vinfp [3]    Outgoing relative velocity to planet vector [km/s]
%  Rpl [1]      Radius of the planet [km]
%  hatm_f [1]   Height of atmosphere [km]
%  input        Input structure
% 
% OUTPUT:
%  dvp [1]      Requested DeltaV at pericenter [km/s] 
%  delta [1]    Turn angle [rad]
%  rp [1]       Radius of pericenter [km]
%  em [1]       Eccentricity of the incoming hyperbola 
%  ep [1]       Eccentricity of the outgoing hyperbola
%  am [1]       Semi-major axis of the incoming hyperbola [km]
%  ap [1]       Semi-major axis of the outgoing hyperbola [km]
%  vpm [1]      Velocity at pericenter of the incoming hyperbola [km/s]
%  vpp [1]      Velocity at pericenter of the outgoing hyperbola [km/s]
%  deltam [1]   Turn angle of the incoming hyperbola [rad]
%  deltap [1]   Turn angle of the outgoing hyperbola [rad]
% 
% CONTRIBUTORS:
%  Lyle Campbell
%  Giulio Pacifici 
%  Luca Rizzieri 
%  Davide Sisana
% 
% VERSIONS:
%  24-11-2019: First version
%

delta = acos(dot(vinfm,vinfp)/norm(vinfm)/norm(vinfp));

emp = @(rp,vinf) 1 + rp*vinf^2/input.k_f;
deltamp = @(rp,vinf) 2*asin(1/emp(rp,vinf));

f = @(rp) delta - 0.5*deltamp(rp,norm(vinfm)) ...
    -0.5*deltamp(rp,norm(vinfp));

if (f(Rpl + hatm_f)<=0)&&(f(1e10)>0)
    
    options = optimset('TolX',1e-14);
    
    rp = fzero(f,[Rpl + hatm_f, 1e10],options);
    
    vpp = sqrt(norm(vinfp)^2 + 2*input.k_f/rp);
    vpm = sqrt(norm(vinfm)^2 + 2*input.k_f/rp);
    
    dvp = norm(vpp - vpm);

else
    dvp = NaN;
    rp = NaN;
    
end

em = emp(rp,norm(vinfm));
ep = emp(rp,norm(vinfp));

deltam = deltamp(rp,norm(vinfm));
deltap = deltamp(rp,norm(vinfp));

am = rp/(1-em);
ap = rp/(1-ep);
end


