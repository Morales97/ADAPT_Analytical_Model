%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate state probability (number of paquets in queue)
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fileID = fopen('state_result_3way_50n_200us_3.txt');
formatSpec = '%i %i'; 
dims = [2 Inf];
data = fscanf(fileID, formatSpec, dims);
data = data';

max_state = 5;
state_prob = zeros(max_state, 1);
for h = 1:max_state
    state_prob(h) = length(data(data(:,2) == (h-1)));
    state_prob(h) = state_prob(h) / length(data);
end
state_prob(max_state) = state_prob(max_state) + (length(data(data(:,2) > max_state)) / length(data));

state_prob
            
            
            
            
