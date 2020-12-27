%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute average Tcycle of a simulation
% Plot a histogram showing Tcycle distribution
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fileID = fopen('BS_result_3way_50n_500us_1.txt');
formatSpec = '%i'; 
dims = [1 Inf];
data = fscanf(fileID, formatSpec, dims);
data = data';

cycles = diff(data);
%length(cycles)
Tcycle = mean(cycles)

figure()
histogram(cycles, 100)
title('Tcycle distribution')
xlabel('Tcycle [ns]')
hold on
xline(Tcycle, 'red', 'mean', 'LineWidth', 1.5)

figure()
histogram(cycles*1e-3, 100)
title('Tcycle distribution')
xlabel('Tcycle [\mus]')
hold on
xline(Tcycle*1e-3, 'red', 'mean', 'LineWidth', 1.5)


set(gca, 'FontSize', 14, 'LineWidth', 1)

