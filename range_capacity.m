%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the analytical range for every point [node_density, tia]
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Unused - Find stability limit for each scenario
% % 15 m, 8PSK
% radius = 15;
% Nsec = 30;
% data_rate = 157e9;
% 
% Tia_arr = [20, 30, 40, 50, 60, 70, 80, 100, 125];
% %Tia_arr = [20, 30, 40, 50];
% nodesInit = 10;   
% 
% analytical_stability_limit_simplified;
% % This is computing the point where avgT > Tia. It is not exactly the same
% % as the point where rho > 1, there will be a little difference compared to
% % theoretical approach
% plot(nodeDensities, Tia_arr, 'x', 'LineWidth', 3, 'DisplayName', '15m analytical')
% fit = polyfit(nodeDensities, Tia_arr, 2);
% plot(x_values, polyval(fit, x_values),'DisplayName', 'line fit')
% 
% % Theoretical approach
% [Ttx, Tskip] = calc_times(radius, Nsec, data_rate);
% plot(x_values, polyval([Ttx*pi*radius^2, Nsec*Tskip], x_values),'DisplayName', 'theo')
% legend('show')

%% Ranges comparison
x_values = 0:0.005:1;

radius_arr = [7.5, 10, 13, 15, 18];
%radius_arr = [7.5];
Nsec = 30;
range_mcs = [7.5, 16, 18];
data_rate = [315, 210, 157]*1e9;
%range_mcs = [18];
%data_rate = [157e9];
rho_aimed = 0.11; % 1 for stability limit, 0.1 for high throughput

figure()
for g = 1:length(radius_arr)
    radius = radius_arr(g);
    [Ttx, Tskip] = calc_times(range_mcs, radius, data_rate);
    hold on
    plot(x_values, polyval([Ttx*pi*radius^2, Nsec*Tskip/rho_aimed], x_values), 'DisplayName',[num2str(radius) ' m, theo'])
end
xlabel('Node density [nodes/m^2]');
ylabel('Tia [us]');
ylim([0 600])
legend('show')

%% Range for every point (Tia, node density) for a certain target rho
Nsec = 30;
radius = 25;
data_rate = [315, 210, 157, 105.3]*1e9;
%data_rate = 210e9;
%range_mcs = [9.5, 17.6, 18];
range_mcs = [7.5, 15.4, 18, 34];

Tia = 0:3:1000;
node_density = 0.03:0.01:1;
alpha = 0.5;
range = zeros(length(Tia), length(node_density));

% Target rho
rho = 0.1;

for i = 1:length(Tia)
    for j = 1:length(node_density)
        
        % 1. Calculate using highest order modulation
        [Ttx_noProp, Tskip_noProp] = calc_times_noProp(data_rate(1));
        
        a = rho * node_density(j) * pi * alpha * 3.33e-3;
        b = rho * node_density(j) * pi * Ttx_noProp;
        c = Nsec * 6.66e-3;
        d = Nsec * Tskip_noProp - rho * Tia(i);
        roots_ = roots([a b c d]);
        
        range(i,j) = max(real(roots_));
        if range(i,j) < 0
            range(i,j) = 0;
        end
        
        if range(i,j) > range_mcs(1)    
            % 2. If range is more than highest order modulation, calculate using also 2nd highest MCS
            [Ttx_noProp_2, Tskip_noProp] = calc_times_noProp(data_rate(2));
            b = rho * node_density(j) * pi * Ttx_noProp_2;
            d = Nsec * Tskip_noProp - rho * Tia(i) + rho * node_density(j) * pi * range_mcs(1)^2*(Ttx_noProp - Ttx_noProp_2);
            roots_ = roots([a b c d]);
        
            range(i,j) = max(real(roots_));
        end
        
        if range(i,j) > range_mcs(2)    
            % 3. Same procedure for the 2nd to 3rd MCS boundary
            [Ttx_noProp_3, Tskip_noProp] = calc_times_noProp(data_rate(3));
            b = rho * node_density(j) * pi * Ttx_noProp_3;
            d = Nsec * Tskip_noProp - rho * Tia(i) + rho * node_density(j) * pi * (range_mcs(1)^2*(Ttx_noProp - Ttx_noProp_2) + range_mcs(2)^2 * (Ttx_noProp_2 - Ttx_noProp_3));
            roots_ = roots([a b c d]);
        
            range(i,j) = max(real(roots_));
        end
        
        if range(i,j) > range_mcs(3)    
            % 4. 3nrt to 4th MCS boundary
            [Ttx_noProp_4, Tskip_noProp] = calc_times_noProp(data_rate(4));
            b = rho * node_density(j) * pi * Ttx_noProp_4;
            d = Nsec * Tskip_noProp - rho * Tia(i) + rho * node_density(j) * pi * (range_mcs(1)^2*(Ttx_noProp - Ttx_noProp_2) + range_mcs(2)^2 * (Ttx_noProp_2 - Ttx_noProp_3)  + range_mcs(3)^2 * (Ttx_noProp_3 - Ttx_noProp_4));
            roots_ = roots([a b c d]);
        
            range(i,j) = max(real(roots_));
        end
        
        if range(i,j) > radius
            range(i,j) = radius;
        end
    end
end

figure()
surf(node_density, Tia, range, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2)
ylim([min(Tia) max(Tia)])
xlim([min(node_density) max(node_density)])
colormap(jet)
c = colorbar;
brighten(0.3)
ylabel(c, 'Range [m]', 'FontSize', 12)
xlabel('Node density [nodes/m^2]')
ylabel('Tia [\mus]')
zlabel('Range [m]')
set(gca, 'FontSize', 12)

%title('Maximum range to keep system stable [m]')

%% Rho comparison
x_values = 0:0.005:1;

Nsec = 30;
range_mcs = [7.5, 15.4, 18];
%radius = max(range_mcs);
radius = 18;
data_rate = [315, 210, 157]*1e9;

[Ttx, Tskip] = calc_times(range_mcs, radius, data_rate);

%figure()
hold on
plot(x_values, polyval([Ttx*pi*radius^2, Nsec*Tskip], x_values),'DisplayName','rho=1')
plot(x_values, polyval([Ttx*pi*radius^2, Nsec*Tskip/0.5], x_values),'DisplayName','rho=0.5')
plot(x_values, polyval([Ttx*pi*radius^2, Nsec*Tskip/0.25], x_values),'DisplayName','rho=0.25')
plot(x_values, polyval([Ttx*pi*radius^2, Nsec*Tskip/0.1], x_values),'DisplayName','rho=0.1')
legend('show')
xlabel('Node density [nodes/m^2]');
ylabel('Tia [us]');
ylim([0 600])
%xlim([0.01 0.2])


%% Calculate Ttx and Tskip [us]

function [Ttx, Tskip] = calc_times(range_mcs, radius, data_rate)

Tcts = 20*8 / min(data_rate) * 1e9;    % ns
Trts = Tcts;
Tack = Tcts;
Tprop_max = radius * 3.33;    % ns
Tprop = Tprop_max/2;
Tbo_start = 10;
packet_size = 65000;

% Limit range to radius
range_mcs(range_mcs > radius) = radius;

prob_range = zeros(1, length(range_mcs));
Tdata = 0;
for i = 1:length(prob_range)
    if i == 1
        prob_range(1) = range_mcs(1)^2 / max(range_mcs)^2;
    else
        prob_range(i) = (range_mcs(i)^2 - range_mcs(i-1)^2) / max(range_mcs)^2; 
    end
    Tdata = Tdata + prob_range(i) * (packet_size * 8 / data_rate(i)) * 1e9;
end

% % Node distribution
% beamwidth = 360 / Nsec; % degrees
% A = beamwidth / 360 * pi * radius^2; % [m^2]
% lambda_A = nodes / (pi * radius^2);  % [nodes/m^2]
% 
% Pi_in_A_aux = zeros(1, nodes);
% for i = 0:nodes-1
%     Pi_in_A_aux(i+1) = (lambda_A * A)^i / factorial(i) * exp(-lambda_A * A);
% end
% Pi_in_A = Pi_in_A_aux(1:Nmax_nodes);
% Pi_in_A(Nmax_nodes + 1) = sum(Pi_in_A_aux(Nmax_nodes+1:nodes)); 

Ttx1 = Tcts + Tdata + Tack + Tprop * 2; % first node of the sector
Ttx2 = Tdata + Tbo_start + Tcts + Tack;   % rest of nodes of the sector

%Ttx = (Ttx1 * (Pi_in_A(2)) + Ttx2 * (1- Pi_in_A(1) - Pi_in_A(2))) / (1-Pi_in_A(1));
Ttx = (Ttx1 * 0.5 + Ttx2 * 0.5) * 1e-3; % Assumption: half of the transmissions are the first of the sector and the other half not. Otherwise, Ttx would depend on the number of nodes as well. 

Tskip = (Tprop_max * 2 + Tbo_start + Tcts + Trts + 1) * 1e-3;  

end


%% 
function [Ttx, Tskip] = calc_times_noProp(data_rate)

Tcts = 25*8 / data_rate * 1e9;    % ns
Trts = Tcts;
Tack = Tcts;
Tbo_start = 10;
packet_size = 65000;

Tdata = (packet_size * 8 / data_rate) * 1e9;

Ttx1 = Tcts + Tdata + Tack; % first node of the sector
Ttx2 = Tdata + Tbo_start + Tcts + Tack;   % rest of nodes of the sector

Ttx = (Ttx1 * 0.5 + Ttx2 * 0.5) * 1e-3;
Tskip = (Tbo_start + Tcts + Trts + 1) * 1e-3;  

end