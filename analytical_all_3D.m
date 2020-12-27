%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical results for all Node-Tia scenarios for a certain configuration
% Model follows the Mathematical framework developed
%
% Set the parameters according to the desired configuration
%   - Nsec: number of sectors
%   - range_mcs: Max distance for every MCS used, from higher to lower modulation order
%   - data_rate: Data rate for every MCS used, from higher to lower modulation order
%   - radius: maximum range. By default max(range_mcs), but can be configured
%   to limit the range to a shorter distance
%
% Daniel Morales - June 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------------
% -------- PARAMETERS ---------
% -----------------------------

% CONFIG 20
Nsec = 30;
range_mcs = [7.5, 16, 18];
data_rate = [315, 210, 157]*1e9;
radius = max(range_mcs);

% CONFIG 29
% Nsec = 30;
% range_mcs = [7.5];
% data_rate = [315]*1e9;
% radius = max(range_mcs);

% Tia_arr = 0:1:300;
% node_density = 0.001:0.002:0.2;
% nodes_arr = node_density * pi * radius^2;
Tia_arr = 100:5:300;
% nodes_arr = 1:1:200;
nodes_arr = 25:5:200;
node_density = nodes_arr / (pi * radius^2);

%% Estimate rho (2D colormap)
[Ttx, Tskip] = calc_times(range_mcs, radius, data_rate);
rho = zeros(length(Tia_arr), length(nodes_arr));

for i = 1:length(Tia_arr)
    for j = 1:length(nodes_arr)
        rho(i,j) = Nsec * Tskip / (Tia_arr(i) - nodes_arr(j) * Ttx);
        if rho(i,j) > 1 || rho(i,j) < 0
            rho(i,j) = 1;
        end
    end
end

%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.7 0.9]); 
%subplot(1,2,1)
figure()
surf(nodes_arr / (pi * radius^2), Tia_arr, rho, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(2)
ylim([min(Tia_arr) max(Tia_arr)])
xlim([min(node_density), max(node_density)])
colormap(jet)
colorbar
xlabel('Node density [nodes/m^2]')
ylabel('Tia [us]')
title('rho')


%% ESTIMATE THROUGHPUT
% WARNING. May take long to complete
% Double check parameters in 'analytical_single.m'

%Tia_arr = 50:10:600;
%nodes_arr = 20:2:150;
node_density = nodes_arr / (pi * radius^2);

Throughput = zeros(length(Tia_arr), length(nodes_arr));
Tcycle_matrix = zeros(length(Tia_arr), length(nodes_arr));


for q = 1:length(Tia_arr)
    Tia = Tia_arr(q);
    for r = 1:length(nodes_arr)
        nodes = nodes_arr(r);
        analytical_single
        Throughput(q,r) = S * 1e9;
        %Tcycle_matrix(q,r) = Tcycle * 1e-3;
        %Tcycle = 0;
        S = 0;
   end
end

%% Throughput with PDF, 2D plot

figure()
%subplot(2,2,3)
surf(nodes_arr / (pi * radius^2), Tia_arr, Throughput*1e-9, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
view(2)
ylim([min(Tia_arr) max(Tia_arr)])
xlim([min(node_density), max(node_density)])
colormap(jet)
c = colorbar;
caxis([0 200])
xlabel('Node density [nodes/m^2]')
ylabel('Tia [us]')
ylabel(c, 'Throughput [Gbps]')
%title('Throughput estimation')
grid off
set(gca, 'FontSize', 15, 'LineWidth', 1)

%% Throughput with PDF, 3D plot
figure()
%subplot(2,2,4)
surf(nodes_arr / (pi * radius^2), Tia_arr, Throughput, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
view(3)
ylim([min(Tia_arr) max(Tia_arr)])
xlim([min(node_density), max(node_density)])
zlim([0 max([10e10 max(Throughput)])])
colormap(jet)
colorbar
caxis([0 18e10])
xlabel('Node density [nodes/m^2]')
ylabel('Tia [us]')
title('Throughput estimation')

%% Tcycle

figure()
surf(nodes_arr / (pi * radius^2), Tia_arr, Tcycle_matrix, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
view(2)
ylim([min(Tia_arr) max(Tia_arr)])
xlim([min(node_density), max(node_density)])
colormap(jet)
c = colorbar;
caxis([0 10])
xlabel('Node density [nodes/m^2]')
ylabel('Tia [us]')
ylabel(c, 'Tcycle [\mus]')
%title('Throughput estimation')
grid off
set(gca, 'FontSize', 15, 'LineWidth', 1)


%% Function used to calculate the transmission and skip time of each config
% Calculate Ttx and Tskip [us]

function [Ttx, Tskip] = calc_times(range_mcs, radius, data_rate)

Tcts = 25*8 / min(data_rate) * 1e9;    % ns
Trts = Tcts;
Tack = Tcts;
Tprop_max = radius * 3.33;    % ns
Tprop = Tprop_max/2;
Tbo_start = 10;

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
    Tdata = Tdata + prob_range(i) * (15000 * 8 / data_rate(i)) * 1e9;
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