%% Parameters
% Select parameters of files to merge
% Ex: node_num = [10, 25]; time_of_enqueue = [22]; will merge all files
% with 10 or 25 nodes and 22us of enqueue time

num_nodes = [50];
tia = [200, 300, 400, 500, 600, 800, 1000];
%num_nodes = [10, 14, 18, 22, 26];
%num_nodes = [40, 55, 70, 85, 100];
%num_nodes = [50, 100, 150, 200];
%num_nodes = [30, 60, 90, 120, 150];
%num_nodes = [100, 150, 200];
%num_nodes = [40];
%tia = [200, 100, 75, 50];     % Microseconds
%tia = [50, 75, 100, 200];
%tia = [100, 200, 300, 500];
%tia = [300, 500];
%radius = 21;          

%% Check folder
path = './';
% check to see if called the correct folder 
if exist(path, 'dir')~= 7
   Message = sprintf('Error: The following folder does not exist:\n%s', path);
   uiwait(warndlg(Message));
   return;
end

%% 
thMatrix = zeros(length(num_nodes), length(tia));
timeMatrix = zeros(length(num_nodes), length(tia));
PdisMatrix = zeros(length(num_nodes), length(tia));
avgThMatrix = zeros(length(num_nodes), length(tia));
avgTimeMatrix = zeros(length(num_nodes), length(tia));
avgPdisMatrix = zeros(length(num_nodes), length(tia));
%%
for i = 1 : length(num_nodes)
    for j = 1 : length(tia)
        sname = sprintf('result_0way_%un_%uus_1.txt', num_nodes(i), tia(j));
        filePattern = fullfile(path, sname);
        FileList = dir(filePattern);

        res = zeros(1,length(FileList));
        for k = 1:length(FileList)
            baseFileName = FileList(k).name;
            
            % Compute metrics for this scenario
            fileID = fopen(baseFileName);
            formatSpec = '%i %i %f %i %i'; 
            dims = [5 Inf];
            data = fscanf(fileID, formatSpec, dims);
            %fclose(fullfile(path, baseFileName));
            [thr, time, time_perm, Pdis] = computeMetrics(data,num_nodes(i));
            thMatrix(i,j,k) = thr;
            timeMatrix(i,j,k) = time_perm;
            PdisMatrix(i,j,k) = Pdis;
        end
        sum(thMatrix(i,j,:));
        avgThMatrix(i,j) = sum(thMatrix(i,j,:)) / length(FileList);
        avgTimeMatrix(i,j) = sum(timeMatrix(i,j,:)) / length(FileList);
        avgPdisMatrix(i,j) = sum(PdisMatrix(i,j,:)) / length(FileList);
    end
end

% figure()
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.7 0.9]);
% 
% % Average Throughput vs tia
% subplot(2,2,1)
%subplot(2,1,2)
for i=1:length(num_nodes)
    hold on
    plot(tia, avgThMatrix(i,:), '-o', 'DisplayName', [num2str(num_nodes(i)) ' nodes sim'])
    title(['Avg Throughput'])
    xlabel('Inter arrival time [us]')
    %ylim([0, 4e10])
end
legend('show')

% % Average Throughput vs node desnity
% subplot(2,2,2)
% for i=1:length(tia)
%     hold on
%     plot(num_nodes ./ (pi * radius ^2), avgThMatrix(:,i), '-o', 'DisplayName', [num2str(tia(i)) ' us'])
%     title(['Avg Throughput'])
%     xlabel('Node density [nodes/m^2]')
%     ylim([0, 4e10])
% end
% legend('show')
 
% Average Packet time
% limit = length(tia);
% %subplot(2,2,3)
% subplot(2,2,1)
% for i=1:length(num_nodes)
%     plot(tia(1:limit), avgTimeMatrix(i,1:limit), '-o', 'DisplayName', [num2str(num_nodes(i)) ' nodes sim'])
%     hold on
%     title(['Avg Packet Time'])
%     xlabel('Inter arrival time [us]')
% end
% %plot(tia, tia, 'DisplayName', ['Stability condition'])
% legend('show')
% 
% % Average Packet time
% subplot(2,2,4)
% for i=1:length(tia)
%     plot(num_nodes ./ (pi * radius ^2), avgTimeMatrix(:,i), '-o', 'DisplayName', [num2str(tia(i)) ' us'])
%     hold on
%     title(['Avg Packet Time'])
%     xlabel('Node density [nodes/m^2]')
% end
% yline(max(tia),'-.', [num2str(max(tia)) ' us'], 'LabelHorizontalAlignment', 'left');
% legend('show')
% legend('Location', 'northwest')


%% Compute metrics 
% data: Matrix of results for the node. 
% node: Number of nodes

function [throughput, time, time_perm, discard_rate] = computeMetrics(data, num_nodes)

data = data';
packet_size = data(1,2);
throughput = 0;
time = 0;

for j = 1:num_nodes
    th_aux = 0;
    time_aux = 0;
    node_data = data(data(:,1)==j,:);           % Select only data from node j
    succ_data = node_data(node_data(:,4)==1,:); % Select successfull packets 
    num_succ = length(succ_data(:,1));
    for i = 1:num_succ
        th = packet_size * 8 / (succ_data(i,3) * 1e-9);
        th_aux = th_aux + th;
        time_aux = time_aux + (succ_data(i,3) * 1e-3);
    end
    th_node = th_aux / num_succ;
    time = time + time_aux / num_succ;
    if num_succ == 0
        th_node = 0;
    end
    throughput = throughput + th_node;
end

% average time on permanent phase
succ = data(data(:,4)==1, :);
time_perm = mean(succ(round(length(succ)/10) : length(succ), 3)) * 1e-3;

throughput = throughput / num_nodes;
time = time / num_nodes;
discard_rate = length(data(data(:,5)==1,:)) / length(data);
end