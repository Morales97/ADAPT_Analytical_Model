%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate analytical results for a set of scenarios

% ATTENTION
% Make sure to check that parameters in 'analytical_single.m' are correct
% Specially 'range_mcs', which has to indicate the range of each modulation
% Also, comment 'num_nodes' and 'tia' in 'sim_results' and 'sim_tcycle_multi'

% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------------
% -------- PARAMETERS ---------
% -----------------------------

num_nodes = [25, 50];
%num_nodes = [100]
%tia = [200, 100, 50, 25];
tia = [1000, 800, 600, 400, 300, 200, 175, 150];%, 125, 100]; % Microseconds

radius =  18;

%% -------------------------------------------
% --- Run analytical model for all configs ---
% --------------------------------------------

avgTimeTheo = zeros(length(num_nodes), length(tia));
avgThTheo = zeros(length(num_nodes), length(tia));
avgTcycleTheo = zeros(length(num_nodes), length(tia));

for k = 1 : length(num_nodes)
    nodes = num_nodes(k);
    
    for j = 1 : length(tia)
        Tia = tia(j);
 
        analytical_single
        avgTimeTheo(k, j) = avgT;
        avgThTheo(k, j) = S*1e9;
        avgTcycleTheo(k, j) = Tcycle;
    end
    
    %figure()
    subplot(2,1,1)
    hold on
    plot(tia, avgTcycleTheo(k,:).*1e-3, '-o', 'DisplayName', [num2str(num_nodes(k)) ' nodes analytical']);
    title('Avg Tcycle')
    xlabel('Inter arrival time [us]')
    ylabel('Tcycle [us]')
    
    subplot(2,1,2)
    hold on
    plot(tia, avgThTheo(k,:), '-o', 'DisplayName', [num2str(num_nodes(k)) ' nodes analytical']);
    title('Avg Throughput')
    xlabel('Inter arrival time [us]')
    %ylim([0 10e10])
    ylabel('bps')
    legend('Location', 'southeast')
    
end
legend('show')

%% Compare analytical model to ns3 simulations - need results files in the correct folders
%sim_tcycle_multi
%sim_results