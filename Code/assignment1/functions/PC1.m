function  PC1(t_d,t_f,x_ga,x_gs, input)
% PC1 pork chop plot matrices for the first transfer
%
% Function to compute the pork chop plot for the first transfer.
% 
% PROTOTYPE:
%  PC1(t_d,t_f,x_ga,x_gs,input)
% 
% INPUT:
%  t_d          Departure time array
%  t_f          Flyby time array
%  x_ga         Point of minimum from ga
%  x_gs         Point of minimum from grid search
%  input        Input structure
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
%  08-12-2019: First version
% 

% Load constants:

% Given Planets
id_d = input.id_d;                                   % Neptune
id_f = input.id_f;                                   % Mars

ksun = input.ksun;               % Sun gravity constant

% Make grids:

[TD1,TF1] = meshgrid(t_d,t_f);
DV1 = TD1+NaN;

% For loops:

for j = 1:length(t_d)
    
    for i = 1:length(t_f)
        
        % Check1: t_d<t_f<t_a must be True
        if  (t_d(j)<=t_f(i))
            kep_d = uplanet(t_d(j), id_d);
            kep_f = uplanet(t_f(i), id_f);
            
            [rr_d,vv_d] = par2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
            [rr_f,vv_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
            
            [~,~,~,err1,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f(i)-t_d(j))*24*3600,ksun,0,0,0,0);
            
            % Check2: there must be no errors in lambertMR
            if (err1==0)
                DV1(i,j) = norm(vv_d - vt1_i');
            end
            
        end
    end
end

% Plots:

time_shift = 730486.5; % Time shift for plots

figure,
grid on
hold on
contour(TD1+time_shift,TF1+time_shift,DV1,[0:0.5:20],'HandleVisibility','off');
datetick('x','yyyy mmm dd')
xtickangle(45)
datetick('y','yyyy mmm dd')
ytickangle(45)
hcb = colorbar;
set(get(hcb,'Title'),'String','\DeltaV [km/s]');
title('Neptune - Mars')
xlabel('Departure Time')
ylabel('Flyby Time')

plot(x_ga(1)+time_shift,x_ga(2)+time_shift,'ro','linewidth',2)
plot(x_gs(1)+time_shift,x_gs(2)+time_shift,'ko','linewidth',2)

legend('ga Result','Grid Search Result')
hold off

