%% Upload data

fileID = fopen('result_2way_50n_200us.txt', 'r');

formatSpec = '%i %i %f %i %i'; 
dims = [5 Inf];

data = fscanf(fileID, formatSpec, dims);
data = data';
numNodes = max(data(:,1));

% Remove transitory
%data = data(round(length(data)/10):length(data), :);

%% Process data
discarded = data(data(:,5)==1,:);
successful = data(data(:,4)==1,:);

P_discard = length(discarded) / length(data);
avg_delay_us = mean(successful(:,3)) * 1e-3;
avg_delay_perm = mean(successful(round(length(successful)/10) : length(successful), 3)) * 1e-3;

%% Graph tx time
avg_tx_time = zeros(length(successful), 1);
avg_tx_time(1) = successful(1,3);
for i=2:length(successful)
    avg_tx_time(i) = (avg_tx_time(i-1)*(i-1) + successful(i,3)) / i;
end

figure()
h = histogram(successful(:,3),300);

% wait_0 = successful((successful(:,3) <= 1867 + 1753),:);
% wait_1 = successful((2*1867 + 1753 >= successful(:,3)),:);
% wait_1 = wait_1((wait_1(:,3) > 1867 + 1753),:);
% wait_2 = successful(3*1867 + 1753 >= successful(:,3),:);
% wait_2 = wait_2((wait_2(:,3) > 2*1867 + 1753),:);
% 
% tface = successful(:,3) - 1753;
% wait_0 = sum(tface < 967)
% wait_1 = sum(tface >= 967)
% max_tface = max(tface)
% min_tface = min(tface)
% avg_tface = mean(tface)
% figure()
% histogram(tface, 300)
%xlim([0 2000])

% figure()
% p_time = h.Values ./ length(h.Data);
% th_using_pdf = sum(15000 * 8 ./(h.BinEdges(1:length(h.BinEdges)-1) + h.BinWidth) .* p_time)
% title('Packet time')
% xlim([0 50000])

figure()
subplot(2,2,1)
plot(1:length(successful), successful(:,3)*1e-9)
title('Total time for the Nth packet')
% title('25 nodes - 50 us')
%sgtitle('1-way. config 29. 10n 200us seed1')

hold on
plot(1:length(successful), avg_tx_time*1e-9)


%figure()
subplot(2,2,2)
plot(1:length(successful), avg_tx_time*1e-9)
title('Average Tx time vs packets successful')

discard_accum = zeros(length(data), 1);
discard_accum(1) = data(1,5);
for i=2:length(data)
    discard_accum(i) = discard_accum(i-1) + data(i,5);
end

% figure()
% plot(1:length(data),discard_accum)
%% Processing per node

% this step is computationally inefficient and large simulations are stuck
% here



data_node = {};
succ_node = {};
packet_size = data(1,2);

for j = 1:numNodes
    data_node{j} = data(data(:,1)==j,:);
    succ_node{j} = data_node{j}(data_node{j}(:,4)==1,:);
end


th_per_node = zeros(numNodes,1);
P_disc_per_node = zeros(numNodes,1);
succ_per_node = zeros(numNodes,1);
disc_per_node = zeros(numNodes,1);
for j = 1:numNodes
    th_aux = 0;
    for i = 1:length(succ_node{j}(:,1))
        th = packet_size * 8 / (succ_node{j}(i,3) * 1e-9);
        th_aux = th_aux + th;
    end
    P_disc_per_node(j) = sum(data_node{j}(:,5)==1) / length(data_node{j});
    th_per_node(j) = th_aux / length(succ_node{j}(:,1));
    succ_per_node(j) = length(succ_node{j}(:,1));
    disc_per_node(j) = length(data_node{j}(:,1)) - succ_per_node(j);
end

th_per_node(isnan(th_per_node)) = 0;    % replace NaN for 0
P_disc_per_node(isnan(P_disc_per_node)) = 1;

throughput = sum(th_per_node);
avg_th_per_node = throughput / numNodes
%succ_node = succ_node';
total_packets = sum(succ_per_node);
%% Graphs per node
avg_Ttx_per_node = zeros(numNodes,1);
for j = 1:numNodes
    avg_Ttx_per_node(j) = mean(succ_node{j}(:,3)*1e-3);
end

%figure()
subplot(2,2,3)
for j = 1:numNodes
    plot(1:size(succ_node{j},1), succ_node{j}(:,3)*1e-9)
    hold on
end
title('Time for the Nth packet per node')


subplot(2,2,4)
batData = [succ_per_node disc_per_node];
bar(batData, 'stacked')
title('Packets succ/disc by node')
    
