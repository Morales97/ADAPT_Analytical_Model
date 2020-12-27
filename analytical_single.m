%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate analytical results for a single scenario 
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------------
% -------- ASSUMPTIONS --------
% -----------------------------

% Assumption 1: RTS don't collide
% Assumption 2: Nodes only use the sector where they have more power
% Assumption 3: Tprop is approximated as Tprop_max/2

% BS Logic: ALL_RTS_MCS - Answers to all RTS. In practice, this means
% that in one Tcycle, every node will have one TX opportunity
% (independently of distance to BS and number of nodes in the sector)


%% ----------------------------
% -------- PARAMETERS ---------
% -----------------------------

% Change according to simulation settings
%config 20
range_mcs = [9.5, 17.6, 18];       % Range in m of each MCS used, from higher order MCS to lower order
data_rate = [315, 210, 157]*1e9;   % [bps] For each MCS used, from higher order MCS to lower order
Nsec = 30;

% config 29
% range_mcs = [7.5];        % Range in m of each MCS used, from higher order MCS to lower order
% data_rate = [315]*1e9;    % [bps] For each MCS used, from higher order MCS to lower order
% Nsec = 30;

% Comment nodes and tia if using 'analytical_all'
% nodes = 50;
% Tia = 300;

packet_size = 65000;    % 65 kB of payload
cwMax = 5;
Tslot = 2;              % [ns]
Nmax_sector = 10;       % max nodes per sector, 
Pc = 0;                 % Average collision prob (Pc)
Pc_worst = 0;           % Pc for the worst node

enable_figures = 0;
warning('off')          % Supress warnings to increase execution speed.

%% Autocalculated parameters

% Time
Tcts = 25*8 / min(data_rate) * 1e9;   % ns
Trts = Tcts;
Tack = Tcts;
Tprop_max = max(range_mcs) * 3.33;    % ns
Tprop = Tprop_max * 2 / 3;            % more nodes in longer distances
Tbo_start = cwMax * Tslot;            % ns

prob_range = zeros(1, length(range_mcs));
Tdata_range = zeros(1, length(range_mcs));
Tprop_range = zeros(1, length(range_mcs));
Tdata = 0;           % Average Tdata depends on the MCS used, which depends on the range
for i = 1:length(prob_range)
    if i == 1
        prob_range(1) = range_mcs(1)^2 / max(range_mcs)^2;
    else
        prob_range(i) = (range_mcs(i)^2 - range_mcs(i-1)^2) / max(range_mcs)^2; 
    end
    Tprop_range(i) = 3.33 * range_mcs(i);
    Tdata_range(i) = (packet_size * 8 / data_rate(i)) * 1e9;
    Tdata = Tdata + prob_range(i) * Tdata_range(i);
end

Tskip = Tprop_max * 2 + Tbo_start + Tcts + Trts + 1;  
Ttx1 = Tcts + Tdata + Tack + Tprop * 2;   % first node of the sector
Ttx2 = Tdata + Tbo_start + Tcts + Tack;   % rest of nodes of the sector

% Node distribution
beamwidth = 360 / Nsec; % degrees
A = beamwidth / 360 * pi * max(range_mcs)^2; % [m^2]
lambda_A = nodes / (pi * max(range_mcs)^2);  % [nodes/m^2]

Pi_in_A_aux = zeros(1, nodes);
for i = 0:nodes-1
    Pi_in_A_aux(i+1) = (lambda_A * A)^i / factorial(i) * exp(-lambda_A * A);
end
Pi_in_A = Pi_in_A_aux(1:Nmax_sector);
Pi_in_A(Nmax_sector + 1) = sum(Pi_in_A_aux(Nmax_sector+1:nodes));


%% rho 
syms rho real
Pskip = 0;
for i = 1:Nmax_sector
    Pskip = Pskip + Pi_in_A(i) * (1-rho)^(i-1);
end

% Average Transmission time
Ttx = (Ttx1 * Pi_in_A(2) + Ttx2 * (1- Pi_in_A(1) - Pi_in_A(2))) / (1-Pi_in_A(1)); % It should actually depend on rho too, but the impact is minimal, fair assumption.

% Cycle time
Tcycle = Nsec * Tskip + Ttx * rho * nodes * (1-Pc);

retries = 0;
retries_worst = 0;
n = 2;
for i=1:n   % assuming max n retries
    retries = retries + i * Pc * (1 - Pc)^i;
    retries_worst = retries_worst + i * Pc_worst * (1 - Pc_worst)^i;
end

Tsucc = Tcycle + Tcycle * retries;
Tsucc_worst = Tcycle + Tcycle * retries_worst;

% [-Inf Inf] makes solution real
rho = Nsec*Tskip / (Tia*1e3 - nodes*Ttx);
rho_worst = subs(Tsucc_worst / (Tia*1e3), rho);

Pskip = subs(Pskip, rho);
Pskip = double(Pskip);
Tcycle = subs(Tcycle, rho);
Tcycle = double(Tcycle);
Tsucc = subs(Tsucc, rho);
Tsucc = double(Tsucc);


%% ----------------------------
% ----- ANALYTICAL MODEL ------
% -----------------------------

Tw0 = Tsucc/2;
N = rho / (1-rho)^2;
Tw = N * Tsucc;
Ts = Ttx + Tprop;

% average delay
avgT = double(Tw + Tw0 + Ts);

%% State prob
%sim_states
avg_PstateTheo = [(1-rho), rho * (1-rho), rho^2 * (1-rho), rho^3 * (1-rho)];

%% pdf of Tcycle

if avgT > Tia * 1e3 || avgT < 0 
    S = 0;
elseif rho < 0
    warning('on')
    warning('rho < 0')
else 
    
    Tcycle_mu = nodes * rho * Ttx + Nsec * Tskip;
    Tcycle_sigma = sqrt(nodes * rho * (1-rho) * Ttx);
    limit = round(5 * Tcycle_mu);
    x = -limit:1:limit;
    null_limit = -limit:1:(Nsec * Tskip);
    null_limit2 = 0:1:(Nsec * Tskip);
    pdf_tcycle = normpdf(x, Tcycle_mu, Tcycle_sigma^2);
    cdf_tcycle = normcdf(x, Tcycle_mu, Tcycle_sigma^2);
    null_area = cdf_tcycle(length(null_limit));
    norm_factor = 1/(1-null_area);

    pdf_tcycle_norm = [zeros(1, length(null_limit2)) pdf_tcycle((length(null_limit)+1):length(pdf_tcycle)) * norm_factor];
    if enable_figures
        figure()
        plot(1:1:length(pdf_tcycle_norm),pdf_tcycle_norm)
        title('Estimated PDF of Tcycle')
        xlabel('ns')
        xline(round(Tcycle_mu), 'red', 'Tcycle')
        
        %xlim([0 3000])
    end

    
%% Discrete samples of the probability of Tcycle
    % Compute Tcycle for every number of Tx in the cycle
    % Tcycle_discrete(1): chance of N Tx in the cycle
    % Tcycle_discrete(2): Tcycle for N Tx
       
    Tcycle_discrete = zeros(2,nodes+1);

    for i = 1:nodes+1
        Tcycle_discrete(1,i) = nchoosek(nodes, i-1) * rho^(i-1) * (1-rho)^(nodes-(i-1));
        Tcycle_discrete(2,i) = (i-1) * Ttx + Nsec * Tskip;
        if Tcycle_discrete(1,i) < 0.01 && i > nodes/5
            break   % To avoid excessive computation
        end
    end
    
    if enable_figures
        figure()
        bar(Tcycle_discrete(1,:))
        title('P(n tx in a cycle) estimated using analytical rho')
    end
    
    Tcycle_dist = sum(Tcycle_discrete(1,:) .* Tcycle_discrete(2,:)); % This should match Tcycle 
 
    samples = 5;
    Ts_discrete = zeros(2,samples);
    for i = 1:samples
        Ts_discrete(1,i) = 1/samples;
        Ts_discrete(2,i) = Tskip + Ttx + 2*Tprop * (i-1)/(samples-1);
    end

%% PDF of Tw0 
% Sum of U(0,Tcycle(n)) for every n
    
    pdf_tw0 = zeros(1, limit);

    for i = 1:nodes
        u_limit = min(round(Tcycle_discrete(2,i)), limit);
        uniform_dist = ones(1, u_limit);
        pdf_tw0 = pdf_tw0 + [uniform_dist .* Tcycle_discrete(1,i) zeros(1, length(pdf_tw0) - length(uniform_dist))];
    end

    pdf_tw0_norm = pdf_tw0 ./ sum(pdf_tw0);
        
    if enable_figures
        figure()
        plot(1:1:length(pdf_tw0_norm), pdf_tw0_norm)
        title('PDF of Tw0 (Sum of U(0,Tcycle(n)) for every n)')
        xlabel('ns')
    end
        
%% PDF of Tsector

    % taking into account Adaptive MCS
    pdf_Tsector = zeros(1, round(Tdata*3));
    for i = 1:length(prob_range)
        if i == 1
            pulse = [zeros(1, round(Tcts + Tdata_range(i) + Tack)) 1:1:3*round(Tprop_range(i))];
        else
            pulse = [zeros(1, round(Tcts + Tdata_range(i) + Tack + Tprop_range(i-1))) 3*round(Tprop_range(i-1):1:Tprop_range(i))];
        end
        pulse = pulse ./sum(pulse);
        pdf_Tsector = pdf_Tsector + [prob_range(i) .* pulse zeros(1, length(pdf_Tsector) - length(pulse))];
    end
        
    % probability n=1
    p1 = 0;
    for i = 1:Nmax_sector
       p1 = p1 + Pi_in_A(i+1) * nchoosek(i,1) * rho * (1-rho)^(i-1); 
    end
    
    % probability n=2
    p2 = 0;
    for i = 2:Nmax_sector
       p2 = p2 + Pi_in_A(i+1) * nchoosek(i,2) * rho^2 * (1-rho)^(i-2); 
    end
    
    pdf_Tsector_n2 = conv(pdf_Tsector,pdf_Tsector);
    pdf_Tsector_n2 = pdf_Tsector_n2(1:length(pdf_Tsector));
    
    pdf_n1 = pdf_Tsector ./sum(pdf_Tsector);
    pdf_n2 = pdf_Tsector_n2 ./sum(pdf_Tsector_n2);
    pdf_Tsector = p1 .* pdf_n1 + p2 .* pdf_n2;
    pdf_Tsector = [zeros(1, round(Tskip)) pdf_Tsector];
    pdf_Tsector = pdf_Tsector ./sum(pdf_Tsector);   % normalize
    
    if enable_figures    
        figure()
        plot(1:1:length(pdf_Tsector)+100, [pdf_Tsector zeros(1,100)])
        title('PDF of Tsector assuming n in (1,2)')
        xlabel('ns')
    end
%% Convolute PDFs of Tw0 and Tsector

    pdf_t_conv = conv(pdf_tw0_norm, pdf_Tsector);
    pdf_t_conv = pdf_t_conv ./ sum(pdf_t_conv);
    
    if enable_figures
%         figure()
%         plot(1:1:length(pdf_t_conv), pdf_t_conv)
%         title('Convolution of Tw0 and Tsector')
%         xlabel('ns')
    end
%% PDF of Tw   
    
    p1 = rho * (1-rho);     % Probability Nw = 1. 
    p2 = rho^2 * (1-rho);   % Probability Nw =2
    
%     % PDF of Tcycle based on n and P(n)
%     pdf_Tcyc_tmp = zeros(1, limit);
%     for i = 2:min(length(Tcycle_discrete), 5)
%         value = Tcycle_discrete(2,i);
%         prob = Tcycle_discrete(1,i);
%         term = [zeros(1, round(value)-1) prob zeros(1, length(pdf_Tcyc_tmp) - round(value))];
%         pdf_Tcyc_tmp = pdf_Tcyc_tmp + term;
%     end
%     pdf_Tcyc_tmp = [(1-rho) p1 .* pdf_Tcyc_tmp];
%     
%     %figure()
%     %plot(1:1:length(pdf_Tcyc_tmp), pdf_Tcyc_tmp)
%     pdf_t_conv2 = conv(pdf_t_conv, pdf_Tcyc_tmp);
%     pdf_t_conv2 = pdf_t_conv2 ./ sum(pdf_t_conv2);
%     
%     if enable_figures
%         figure()
%         plot(1:1:length(pdf_t_conv2), pdf_t_conv2)
%         title('Estimated PDF of T2')
%         xlabel('ns')
%         xlim([0 50000])
%     end
     
     
    % PDF of waiting 2 Tcycles
    pdf_2tcyc = conv(pdf_tcycle_norm, pdf_tcycle_norm);
    pdf_2tcyc_norm = pdf_2tcyc ./ sum(pdf_2tcyc); 
    pdf_2tcyc_norm = pdf_2tcyc_norm(1:length(pdf_tcycle_norm));
    
    if enable_figures    
        figure()
        plot(1:1:length(pdf_2tcyc_norm), pdf_2tcyc_norm)
        title('PDF of 2·Tcycle')
        xlabel('ns')
    end
    %pdf_Tw = rho * (1-rho) .* pdf_tcycle_norm;  % only for n_w = 1
    pdf_Tw = p1 .* pdf_tcycle_norm + p2 .* pdf_2tcyc_norm; % n_w = 1,2
    
    pdf_Tw = pdf_Tw(1:length(pdf_tcycle_norm)); 
    % Sum the probability of not waiting (Nw=0)
    pdf_Tw = [(1-rho) zeros(1,length(pdf_tcycle_norm)-1)] + pdf_Tw;
    pdf_Tw_norm = pdf_Tw ./ sum(pdf_Tw);
  
    if enable_figures    
        figure()
        plot(1:1:length(pdf_Tw_norm), pdf_Tw_norm)
        title('PDF of Tw')
        xlabel('ns')
    end
%% Convolute de three pdfs to obtain pdf of T
    pdf_t_conv = conv(pdf_t_conv, pdf_Tw_norm);
    pdf_t_conv = pdf_t_conv ./ sum(pdf_t_conv);
    
    if enable_figures
        figure()
        plot(1:1:length(pdf_t_conv), pdf_t_conv*1.1e5, 'Color', 'red', 'LineWidth', 1.5)
        title('Estimated PDF of T')
        xlabel('ns')
        xlim([0 17000])
    end
    
    S = 0;
    for i=1:1:length(pdf_t_conv)
        S = S + packet_size * 8 / i * pdf_t_conv(i);
    end
    S;    
    
end

