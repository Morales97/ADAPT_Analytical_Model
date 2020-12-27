%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute average Tcycle of a set of simulations
% Plot average Tcycle
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
num_nodes = [50];
tia = [150, 200, 300, 400, 500, 600];     % Microseconds

%% Check folder
path = './other_files/';
% check to see if called the correct folder 
if exist(path, 'dir')~= 7
   Message = sprintf('Error: The following folder does not exist:\n%s', path);
   uiwait(warndlg(Message));
   return;
end

%%
avgTcycleMatrix = zeros(length(num_nodes), length(tia)');

for i = 1 : length(num_nodes)
    for j = 1 : length(tia)
        sname = sprintf('BS_result_3way_%un_%uus_*.txt', num_nodes(i), tia(j));
        filePattern = fullfile(path, sname);
        FileList = dir(filePattern);

        Tcycle = 0;
        for k = 1:length(FileList)
            baseFileName = FileList(k).name;

            fileID = fopen(['other_files/' baseFileName]);
            formatSpec = '%i'; 
            dims = [1 Inf];
            data = fscanf(fileID, formatSpec, dims);
            data = data';

            cycles = diff(data);
            %length(cycles)
            Tcycle = Tcycle + mean(cycles);
        end
            avgTcycleMatrix(i,j) = Tcycle * 1e-3 / length(FileList);    
    end
end

%subplot(2,1,1)
for i=1:length(num_nodes)
    hold on
    plot(tia, avgTcycleMatrix(i,:), '-ob', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', [num2str(num_nodes(i)) ' nodes simulation'])
    xlabel('tia [us]')
    ylabel('Avg Tcycle [us]')
    title('Tcycle (Latency)')
end
legend('show')

%%
yline(3.96, '--b','Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5)
xlim([150 500])
ylim([0 20])
xlabel('Tia [\mus]')
ylabel('Avg Tcycle [\mus]')
title('')
box on
grid on
set(gca, 'LineWidth', 1, 'FontSize', 14)


