function [c,ceq] = rpCon(x, input)
% rpCon constraint for the minimum value of radius of pericenter of the
% flyby manoeuvre.
% 
% Function to enforce the minimum radius of pericenter constraint ( radius
% of Mars + height of Mars atmosphere ) for the minimum search algorithm.
% 
% PROTOTYPE:
%  [c,ceq] = rpCon(x)
% 
% INPUT:
%  x [3,1]        Times array:         x(1) = departure time [MJD2000]
%                                      x(2) = first transfer tof [days]
%                                      x(3) = second transfer tof [days]
% 
% OUTPUT:
%  c              Array of nonlinear inequality ( c < 0 )
%  ceq            Array of nonlinear equality ( ceq = 0 )
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

ceq = 0;

% Initialize times:
x1 = cumsum(x);
t_d = x1(1); t_f = x1(2); t_a = x1(3); 

% Planets id numbers:
id_d = input.id_d;                              % Neptune
id_f = input.id_f;                              % Mars
id_a = input.id_a;                              % Earth

Rpl_f = input.Rpl_f;                            % Mean radius of Mars
hatm_f = input.hatm_f;                          % Height of Mars atmosphere
   
[kep_d,ksun] = uplanet(t_d, id_d);
[rr_d,vv_d] = par2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, id_f);
[rr_f,vv_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,ksun] = uplanet(t_a, id_a);
[rr_a,vv_a] = par2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,err1,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,err2,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

rp = NaN;

% Check2: there must be no errors in lambertMR
if (err1==0)&&(err2==0)
    
    vinfm = vt1_f' - vv_f;
    vinfp = vt2_i' - vv_f;
    
    [~,~,rp] = powGA(vinfm,vinfp,0,0, input);

end

c = Rpl_f + hatm_f - rp; % Inequality variable

% Add constraint tolerance in order to avoid any possible violation of
% constraint:
c = c + input.gaSearch.ConstraintTolerance;
