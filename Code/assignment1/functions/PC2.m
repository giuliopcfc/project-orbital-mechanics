function  PC2(t_f,t_a,x_ga,x_gs, input)
% PC1 pork chop plot matrices for the second transfer
%
% Function to make the pork chop plot for the second transfer.
% 
% PROTOTYPE:
%  PC2(t_f,t_a)
% 
% INPUT:
%  t_f          Flyby time array
%  t_a          Arrival time array
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

% Given Planet
id_a = input.id_a;                                   % Earth
id_f = input.id_f;                                   % Mars

ksun = input.ksun;                                   % Sun gravity constant


% Make time grids:

[TF2,TA2] = meshgrid(t_f,t_a);
DV2 = TF2+NaN;

% For loops:

for j = 1:length(t_f)
    
    for i = 1:length(t_a)
        
        % Check1: t_d<t_f<t_a must be True
        if  (t_f(j)<=t_a(i))
            
            kep_f = uplanet(t_f(j), id_f);
            kep_a = uplanet(t_a(i), id_a);
            
            [rr_f,vv_f] = par2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
            [rr_a,vv_a] = par2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);
            
            [~,~,~,err1,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,(t_a(i)-t_f(j))*24*3600,ksun,0,0,0,0);
            
            % Check2: there must be no errors in lambertMR
            if (err1==0)
                DV2(i,j) = norm(vv_a - vt2_f');
            end
            
        end
    end
end

% Plots:

time_shift = 730486.5; % Time shift for plots

figure,
grid on
hold on
contour(TF2+time_shift,TA2+time_shift,DV2,[0:0.5:20],'HandleVisibility','off');
datetick('x','yyyy mmm dd')
xtickangle(45)
datetick('y','yyyy mmm dd')
ytickangle(45)
hcb = colorbar;
set(get(hcb,'Title'),'String','\DeltaV [km/s]');
title('Mars - Earth')
xlabel('Flyby Time')
ylabel('Arrival Time')

plot(x_ga(2)+time_shift,x_ga(3)+time_shift,'ro','linewidth',2)
plot(x_gs(2)+time_shift,x_gs(3)+time_shift,'ko','linewidth',2)

legend('ga Result','Grid Search Result')
hold off

hold off

