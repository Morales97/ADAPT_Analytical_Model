%% Upload data

fileID = fopen('collisions_result_3way_50n_200us_6.txt', 'r');
%fileID = fopen('combined.txt', 'r');

%formatSpec = '%i %i %i'; 
%size = [3 Inf];
formatSpec = '%i %i'; 
dims = [2 Inf];

data = fscanf(fileID, formatSpec, dims);
data = data';
numNodes = max(data(:,1));
max_retry = max(data(:,2));

% FILL THIS PARAMETER FOR COMPLETENESS
total_succ_packets = 2000;
%total_succ_packets = 9869;
%% Process data

data_node = {};
col_per_node = zeros(1,numNodes)';
for j = 1:numNodes
    data_node{j} = data(data(:,1)==j,:);
    col_per_node(j) = length(data_node{j});
end

num_collisions = length(data)

retry_dist_accum = zeros(max_retry, 1);
retry_dist = zeros(max_retry, 1);
retry_dist_accum(max_retry) = nnz(data(:,2)==max_retry);
retry_dist(max_retry) = retry_dist_accum(max_retry);
for i = 1:(max_retry-1)
    j = max_retry - i;
    retry_dist_accum(j) = nnz(data(:,2)==j);  %nnz: Number of non zeros
    retry_dist(j) = retry_dist_accum(j) - retry_dist_accum(j+1);
end

total_retried_packets = sum(retry_dist);
zero_retries = total_succ_packets - total_retried_packets;
retry_dist_total = zeros(length(retry_dist) + 1, 1);
retry_dist_total = [zero_retries retry_dist'];
avg_retry = sum(retry_dist_total.*(0:max_retry)') / total_succ_packets;

figure()
bar(retry_dist_total, 'XData', 0:length(retry_dist_total)-1)
title('Packet distribution by number of retries')
text(0:length(retry_dist_total)-1,retry_dist_total,num2str(retry_dist_total'),'vert','bottom','horiz','center'); 
box off
    
    