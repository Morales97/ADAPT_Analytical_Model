%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot and calculate range for a certain window
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate range 
% Range to have a 3dB window
% Approximate calcualtion to give a sense of magnitude. Not to be taken
% precisely

TxPower = 20;   % dBm
%CSth = [-48, -45, -42, -38, -32];     % dBm
CSth = [-49, -46, -40.8, -39.2, -33];     % dBm
Gain = [30.59, 28.05, 24.57, 21.04, 18.55] ;   % dBm
%Gain = [52.17, 46.15];
fc = 287e9;
bw = 69e9;
fmax = fc + bw/2;
AbsLoss = 0;    % we assume no absorption loss at this frequency window

%%
m = 1:0.1:100;
RxPower = TxPower + (2.*Gain)' + 10*log10((3e8./(4*pi*fc.*m)).^2) - 3;

plot(m, RxPower)
xlabel('Distance [m]');
ylabel('Rx Power [dBm]');
legend('60 sectors', '45 sectors', '30 sectors', '20 sectors', '15 sectors')
yline(CSth(1),'-.', 'BPSK')
yline(CSth(2),'-.', 'QPSK')
yline(CSth(3),'-.', '8-PSK')
yline(CSth(4),'-.', '16-QAM')
yline(CSth(5),'-.', '64-QAM')

%% Calculate TxPower to reach a certain range (no absorption loss)

MCS_index = 4;  % 1. BPSK; 2. QPSK; 3. 8-PSK; 4. 16-QAM; 5. 64-QAM;
gain_index = 3; % 1. 60 sectors; 2. 45 sectors; 3. 30 sectors; 4. 20 sectors; 5. 15 sectors
range = 15.3;
tx_necessary_power = CSth(MCS_index) - Gain(gain_index)*2 - 10*log10((3e8/(4*pi*fc*range))^2) + 3

%% Calculate Range of each MCS for Adaptive MCS

MCS_index = 4;              % 1. BPSK; 2. QPSK; 3. 8-PSK; 4. 16-QAM; 5. 64-QAM;
gain_index = 3;             % 1. 60 sectors; 2. 45 sectors; 3. 30 sectors; 4. 20 sectors; 5. 15 sectors
beamwidth_3dB = 12;         % According to number of sectors
% phi = 6;                  % angle difference to center of beam [-beamwidth_3dB/2, beamwidth_3dB/2]
Ptx = 20;                   % dBm
max_range = 18;

% gain_dB = -3 * log10(cos(phi/2*pi/180)) / log10(cos(beamwidth_3dB/4*pi/180))
% range = 3e8/(4*pi*fc) / (10^((Ptx - CSth(MCS_index) + gain_dB + Gain(gain_index)*2)/-20))

range = 0;
for phi=0:0.1:beamwidth_3dB/2

    % from TeraSim antenna model
    gain_dB = -3 * log10(cos(phi/2*pi/180)) / log10(cos(beamwidth_3dB/4*pi/180));
    
    range_aux =  3e8/(4*pi*fc) / (10^((Ptx - CSth(MCS_index) + gain_dB + Gain(gain_index)*2)/-20));
    if range_aux > max_range
        range = range + max_range;
    else
        range = range + range_aux;
    end
end

% This is the average range of the MCS from MCS_index. It varies from the
% center of the beam to a side
avg_range = range / length(0:0.1:beamwidth_3dB/2)


