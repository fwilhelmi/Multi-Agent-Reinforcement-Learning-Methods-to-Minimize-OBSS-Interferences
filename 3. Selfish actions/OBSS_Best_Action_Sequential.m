% MASTER THESIS - Minimizing OBSS Interferences using Q-learning and Correlated Equilibria
% Decentralized Q-learning - Thoughfulness Approach
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

clc
clear all
disp('BEST ACTION')
disp('-----------')

% Global variables for the problem
N_WLANs = 10; % Number of coexistent WLANs in the map
NumChannels = 11; % Number of available channels
noise = -100; % Floor noise (dBm)1
B = 20e6; % Available bandwidth (Hz)

% Definition of actions:
%   - 11 possible channels
%   - 2 different levels of CCA
%   - 2 different levels of TxPower
actions_ch = 1:NumChannels;
actions_cca = [-82 -72]; %dBm
actions_tpc = [15 20]; %dBm

R = 500; % Number of repetitions to take the average improvement
maxIter = 1; % Number of iterations for changing network parameters

fairness_start = zeros(1,R); % Fairness at the beginning of the execution at each repetition
fairness_end = zeros(1,R); % Fairness at the end of the execution at each repetition
avg_tpt_start = zeros(1,R); % Avg Throughput at the beginning of the execution at each repetition
avg_tpt_end = zeros(1,R); % Avg Throughput at the end of the execution at each repetition

%CSMA/CA information
CW = 10;
BO=floor(rand()*2^CW);
tMax = 100;

for r=1:R
    disp('-----------')
    txt = ['ROUND ' num2str(r) '/' num2str(R)];
    disp(txt)
    disp('-----------')
    % Generate network of N WLANs and measure Tpt and Fairness
    wlan = GenerateNetwork3D(N_WLANs, NumChannels, B);
    powMat = PowerMatrix(wlan); %
    WLAN_is_transmitting = zeros(1,N_WLANs); % We store here the number of times an AP transmits
    t = 1;
    while t < tMax
        for i=1:N_WLANs, wlan(i).transmitting = 0; end
        BO_values = rand(1,N_WLANs);
        for i=1:N_WLANs
            [a,b]=min(BO_values); % a the value, b its position in the vector
            BO_values(b)=Inf; % to remove it
            wlan(b).transmitting = 1;
            interferences = Interferences(wlan, powMat); % mW
            %BO_values_stas = rand(1,wlan(b).STAs);
            if(pow2db((interferences(b)+db2pow(noise))) < wlan(b).CCA), WLAN_is_transmitting(b)=WLAN_is_transmitting(b)+1;
            else wlan(b).transmitting = 0; end
        end
        t = t+1;
    end
    % CALCULATE SNR AND TPT FOR EACH WLAN
    interferences = Interferences(wlan, powMat); %dBm
    tpt = zeros(1,N_WLANs);
    for i=1:N_WLANs
        wlan(i).snr = powMat(i,i) - pow2db((interferences(i)+db2pow(noise))); % dBm
        wlan(i).tpt = Capacity(wlan(i).BW, db2pow(wlan(i).snr)*(WLAN_is_transmitting(i)/tMax)); % bps
        tpt(i) = wlan(i).tpt;
    end
    avg_tpt_start(r) = mean(tpt);
    fairness_start(r) = JainsFness(tpt);
    txt1 = ['AVERAGE THROUGHPUT (START) = ' num2str(avg_tpt_start(r))];
    disp(txt1)
    txt2 = ['FAIRNESS (START) = ' num2str(fairness_start(r))];
    disp(txt2)
    % Keep track of the tpt on each WLAN for each iteration
    throughput = zeros(N_WLANs, maxIter);
    % Apply the algorithm sequentially for each WLAN
    iteration = 1;  
    while(iteration < maxIter + 1)    
        % Choose best action (channel)
        for i=1:N_WLANs
%             txt = ['WLAN #' num2str(i)];
%             disp(txt)
%             disp('Throughput before applying best action:')
%             disp(wlan(i).tpt)
            [bestCh, Tpt] = bestAction(wlan,i,1,actions_ch);
            % Change parameters according to the action obtained
            wlan(i).channel = bestCh;  
            wlan(i).tpt = Tpt;
            % Repeat for CCA and TPC
            [bestCCA, Tpt] = bestAction(wlan,i,2,actions_cca);
            wlan(i).CCA = bestCCA; 
            wlan(i).tpt = Tpt;
            [bestTPC, Tpt] = bestAction(wlan,i,3,actions_tpc);
            wlan(i).PTdBm = bestTPC;
            wlan(i).tpt = Tpt;
        end
        powMat = PowerMatrix(wlan); %
        WLAN_is_transmitting = zeros(1,N_WLANs); % We store here the number of times an AP transmits
        t = 1;
        while t < tMax
            for i=1:N_WLANs, wlan(i).transmitting = 0; end
            BO_values = rand(1,N_WLANs);
            for i=1:N_WLANs
                [a,b]=min(BO_values); % a the value, b its position in the vector
                BO_values(b)=Inf; % to remove it
                wlan(b).transmitting = 1;
                interferences = Interferences(wlan, powMat); % mW
                %BO_values_stas = rand(1,wlan(b).STAs);
                if(pow2db((interferences(b)+db2pow(noise))) < wlan(b).CCA), WLAN_is_transmitting(b)=WLAN_is_transmitting(b)+1;
                else wlan(b).transmitting = 0; end
            end
            t = t + 1;
        end
        interferences = Interferences(wlan, powMat); %dBm
        for i=1:N_WLANs
            %if wlan(i).transmitting && pow2db(interferences(i)) >= CCA
            wlan(i).snr = powMat(i,i) - pow2db((interferences(i)+db2pow(noise))); % dBm
            wlan(i).tpt = Capacity(wlan(i).BW, db2pow(wlan(i).snr)*(WLAN_is_transmitting(i)/tMax)); % bps
            throughput(i,iteration) = wlan(i).tpt;
%             disp('Throughput after applying best possible action:')
%             disp(wlan(i).tpt)
        end
        iteration = iteration + 1;    
    end
    avg_tpt_end(r) = mean(mean(throughput));
    fairness_end(r) = JainsFness(mean(mean(throughput)));
    txt3 = ['AVERAGE THROUGHPUT (END) = ' num2str(avg_tpt_end(r))];
    disp(txt3)
    txt4 = ['FAIRNESS (END) = ' num2str(fairness_end(r))];
    disp(txt4)
end

disp('-----------')
disp('RESULTS:')
txt5 = ['AVERAGE IMPROVEMENTS IN OVERALL THROUGHPUT = ' num2str(mean(avg_tpt_end)-mean(avg_tpt_start))];
disp(txt5)
txt6 = ['AVERAGE IMPROVEMENTS IN FAIRNESS = ' num2str(mean(fairness_end)-mean(fairness_start))];
disp(txt6)
%l = {}; %legend
%Smooth data with a moving average
smooth_avg_tpt_start = smooth(smooth(avg_tpt_start));
smooth_avg_tpt_end = smooth(smooth(avg_tpt_end));
% Plot the results
plot(1:R,smooth_avg_tpt_start);
hold on;
plot(1:R,ones(1,R)*mean(smooth_avg_tpt_start),'b--');
plot(1:R,smooth_avg_tpt_end);
plot(1:R,ones(1,R)*mean(smooth_avg_tpt_end),'r--');
%     l = [l; ['WLAN ID = ' num2str(i)]];
xlabel('Iteration');
ylabel('Throughput (bps)');
legend('Average Tpt (before)', 'Mean Average Tpt (before)', 'Average Tpt (after)', 'Mean average Tpt (after)');
% axis([1 N 0 max(max(throughput))]);