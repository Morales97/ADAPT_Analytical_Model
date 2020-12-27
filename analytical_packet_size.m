%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot analytical throughput depending on packet size
% Daniel Morales - June 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% CONFIG 20
% Nsec = 30;
range_mcs = [9.5, 17.6, 18];
data_rate = [315, 210, 157]*1e9;
radius = max(range_mcs);

% CONFIG 29
% Nsec = 30;
% range_mcs = [7.5];
% data_rate = [315]*1e9;
% radius = max(range_mcs);


%% Fixed nodes, sweep tia

num_nodes = [50];
%tia = [100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000];
tia = [100, 250, 500, 1000, 2000, 3000, 5000, 10000, 25000];
size = [2000, 15000, 65000, 500000, 1000000, 2000000];


figure()

for k = 1:length(size)
    nodes = num_nodes;
    packet_size = size(k);
    avgThTheo = zeros(1, length(tia));

    for j = 1 : length(tia)
        Tia = tia(j);
        S_2 = 0;
        analytical_single
        avgThTheo(1, j) = S_2*1e9;
    end


    plot(tia*1e-3, avgThTheo(1,:)*1e-9, '-o', 'LineWidth', 1, 'DisplayName', [num2str(packet_size) 'B']);
    hold on
end

grid on
xlabel('Inter arrival time [ms]')
%ylim([0 10e10])
ylabel('Throughput [Gbps]')
legend('show', 'Location', 'southeast')
box on
set(gca, 'FontSize', 15, 'LineWidth', 1)
set(gcf, 'Position', [0 0 550 300])

%% fixed nodes and tia, packet size sweep

nodes = 50;
Tia = 10000; % 10 ms
%size = [15000, 40000, 65000, 150000, 500000, 1000000, 2000000];
size = [10000, 15000, 25000, 40000, 65000, 100000, 200000, 500000, 1000000, 2000000];

% CONFIG 20
Nsec = 30;
range_mcs = [7.5, 16, 18];
data_rate = [315, 210, 157]*1e9;
radius = max(range_mcs);

avgThTheo18 = zeros(1, length(size));
for k = 1:length(size)
    packet_size = size(k);

    analytical_single
    avgThTheo18(1, k) = S_2*1e9;
end


% CONFIG 29
range_mcs = [7.5];
data_rate = [315]*1e9;
radius = max(range_mcs);

avgThTheo7 = zeros(1, length(size));
for k = 1:length(size)
    packet_size = size(k);

    analytical_single
    avgThTheo7(1, k) = S_2*1e9;
end

%% Plot packet size sweep

figure()
%plot(size*1e-3, avgThTheo7(1,:)*1e-9, '-bs', 'LineWidth', 1, 'DisplayName', '7.5m');
semilogx(size*1e-3, avgThTheo7(1,:)*1e-9, '-bs', 'LineWidth', 1, 'DisplayName', '7.5m');
hold on
%plot(size*1e-3, avgThTheo18(1,:)*1e-9, '-s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', '18m');
semilogx(size*1e-3, avgThTheo18(1,:)*1e-9, '-s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', '18m');

grid on
xlabel('Packet size [kB]')
ylabel('Throughput [Gbps]')
legend('show', 'Location', 'southeast')
box on
set(gca, 'FontSize', 15, 'LineWidth', 1)
set(gcf, 'Position', [0 0 550 300])
