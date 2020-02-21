function [dv_min, x] = gridSearch(input)
% gridSearch
% 
% Function to perform the search for minimum deltaV in a grid defined by 
% values of departure time, frist and second transfer arcs time of flight.
% 
% PROTOTYPE:
%  [dv_min, x] = gridSearch(input)
%  
% INPUT:
%  input    Input structure
% 
% OUTPUT:
%  dv_min         Minimum deltaV found in the grid [km/s]
%  x [3,1]        MJD2000 times array for minimum deltaV: x(1) = departure
%                                                         x(2) = flyby
%                                                         x(3) = arrival       
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

% Load constants:

% Time window
ED = input.ED;                                   % Earliest Departure
LA = input.LA;                                   % Latest Arrival

% Planets id numbers:
id_d = input.id_d;                              % Neptune
id_f = input.id_f;                              % Mars
id_a = input.id_a;                              % Earth

% Time arrays:

% TOF limits:
minTOF1 = input.gridSearch.minTOF1; % Minimum time of flight from Neptune to Mars
maxTOF1 = input.gridSearch.maxTOF1; % Maximum time of flight from Neptune to Mars
minTOF2 = input.gridSearch.minTOF2; % Minimum time of flight from Mars to Earth
maxTOF2 = input.gridSearch.maxTOF2; % Maximum time of flight from Mars to Earth

% Time steps:
dt_d = input.gridSearch.dt_d;                 % Departure time
dt_f = input.gridSearch.dt_f;                 % Flyby time
dt_a = input.gridSearch.dt_a;                 % Arrival time

% Time arrays:
t_d = [ED:...
    dt_d:...
    ED+input.TsynME];                            % Departure time
tof1 = [minTOF1:...
    dt_f:...
    maxTOF1];                                    % Flyby time
tof2 = [minTOF2:...
    dt_a:...
    maxTOF2];                                    % Arrival time

% Create results variable:

data_grid = [];

% Compute DeltaV:

% Three nested for cicles
for i = 1:length(t_d)
    
    i/length(t_d)*100                            % Display percentage of progress
    
    % Departure planet position and velocity:
    [kep_d,ksun] = uplanet(t_d(i), id_d);
    [rr_d,vv_d] = par2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
    
    for j = 1:length(tof1)
        
        t_f = t_d(i) + tof1(j); % Flyby time
        
        % Flyby planet position and velocity:
        [kep_f,ksun] = uplanet(t_f, id_f);
        [rr_f,vv_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
        
        [~,~,~,err1,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,tof1(j)*24*3600,ksun,0,0,0,0);
        
        for k = 1:length(tof2)
            
            t_a = t_d(i) + tof1(j) + tof2(k); % Arrival time
            
            % Not feasible times combination
            if tof2(k)+tof1(j)+t_d(i)<LA
                
                % Arrival planet position and velocity:
                [kep_a,ksun] = uplanet(t_a, id_a);
                [rr_a,vv_a] = par2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);
                 
                 [~,~,~,err2,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,tof2(k)*24*3600,ksun,0,0,0,0);
                 
                 dv = NaN; % Initialize dv
                 
                 % Check: there must be no errors in lambertMR
                 if (err1==0)&&(err2==0)
                     
                     % Mars Reference frame
                     vinfm = vt1_f' - vv_f;
                     vinfp = vt2_i' - vv_f;
                     
                     dvp = powGA(vinfm,vinfp,input.Rpl_f, input.hatm_f, input);
                     
                     dv = norm(vv_d - vt1_i')+norm(vv_a - vt2_f') + dvp;
                     
                 end
                
                if not(isnan(dv))         
                    
                    % Save result:  
                    data_grid = [data_grid, [dv; t_d(i);t_f;t_a]];
                          
                end
            end
            
            
        end
    end
end

% Order results wrt deltaV magnitude:
[~,order] = sort(data_grid(1,:));
data_grid = data_grid(:,order);

dv_min = data_grid(1,1);
x = data_grid(2:end,1);
