% MASTER THESIS - Multi-Agent Reinforcement Learning Methods to Minimize OBSS Interferences
% Multi-Agent Multi-Armed Bandits Problem - Upper Confidence Bound
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

clc
clear all

disp('MA-MAB APPROACH - UCB')
disp('-----------------------')
% Global variables for the problem
N_WLANs = 5;
NumChannels = 6; %Number of available channels (from 1 to NumChannels)
noise = -100; % Floor noise (dBm)
B = 20e6; % Available bandwidth (Hz)

N = 1000; % Number of iterations for Q-learning
epsilon = 0.1; % Exploration coefficient
gamma = 0.95; % Discount rate
alpha = 0.1; % Learning rate
% Definition of actions:
%   - 11 possible channels
%   - 2 different levels of CCA
%   - 2 different levels of TxPower
actions_ch = 1:NumChannels;
actions_cca = [-82 -72]; %dBm
actions_tpc = [15 20]; %dBm

R = 20; %Number of iterations for computing the average
% Each state represents an [i,j,k] combination for indexes on "channels", "CCA" and "TxPower"
states_wlan = 1:(size(actions_ch,2)*size(actions_cca,2)*size(actions_tpc,2));

% Variables to keep track of the average throughput and fairness index
% experienced at the start and at the end of the algorithm's executions
avg_tpt_start = zeros(1,R);
fairness_start = zeros(1,R);
avg_tpt_end = zeros(1,R);
fairness_end = zeros(1,R);

% Number of slots to simulate CSMA/CA behavior
time_slots = 100;

% Times we want to inspect the evolution of the tpt and the fairness during Q-learning execution
T = 100;
% Variables to store the evolution of the tpt and the fairness during Q-learning execution
% progressive_tpt = zeros(R,N/T);
% progressive_fairness = zeros(R,N/T);
progressive_tpt = zeros(R,N);
progressive_fairness = zeros(R,N);

for r=1:R % We repeat the experiment during 'R' repetitions for taking the average
    disp('-----------')
    txt = ['ROUND ' num2str(r) '/' num2str(R)];
    disp(txt)
    time_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    disp(time_start)
    disp('-----------')
    % Generate network of N WLANs and measure Tpt and Fairness
    disp('Creating map of OBSS...')
    %wlan = GenerateNetwork3D(N_WLANs,NumChannels,B); 
    load('wlan.mat');
    disp('Channels selected per WLAN')
    for i=1:N_WLANs
        channels(i) = wlan(i).channel;       
    end
    for i=1:NumChannels
        sumChannels(i) = sum(channels==i);
    end
    disp(channels)
    disp('Times a channel is occupied')
    disp(sumChannels)
    powMat = PowerMatrix(wlan); % Compute the power noticed at each WLAN
    count = 1;
    % Define the initial configuration
    for i=1:N_WLANs    
        [~,index_cca] = find(actions_cca==wlan(i).CCA);
        [~,index_tpc] = find(actions_tpc==wlan(i).PTdBm);
        initial_state_indexes = [wlan(i).channel index_cca index_tpc]; %Take indexes of each selected value    
        initial_state(i) = indexes2val(initial_state_indexes(1),initial_state_indexes(2),initial_state_indexes(3),size(actions_ch,2),size(actions_cca,2)); %Convert indexes into a numeric state
    end    
    optimal_reward = Capacity(20e6,db2pow(100)); % Tpt obtained under an excellent SINR and Jain's Fairness Index = 1 
    % Each configuration at each WLAN has a specific reward that is shared with the others
    reward_per_configuration = zeros(N_WLANs,size(states_wlan,2));
    times_played_arm = zeros(N_WLANs,size(states_wlan,2));
    avg_regret_evolution = zeros(1,N);
    regrets = zeros(N_WLANs,N);
    selected_arm = zeros(1,N_WLANs);
    %%%%%%%%%%%%%%%%%%% MA-MAB implementation - START %%%%%%%%%%%%%%%%%%%
    iteration = 1;
    while(iteration <= N)        
        if iteration == 1
            % Compute initial tpt and fairness
            tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
            avg_tpt_start(r) = mean(tpt_aux);
            fairness_start(r) = JainsFness(tpt_aux);
        end
        % Each agent (WLAN) pulls an arm (selects a configuration)
        for i=1:N_WLANs
            % Select best arm - TODO
            selected_arm(i) = ChooseArm(reward_per_configuration(i,:),epsilon);
            times_played_arm(i,selected_arm(i)) = times_played_arm(i,selected_arm(i)) + 1;
            [a b c] = val2indexes(selected_arm(i),NumChannels,size(actions_cca,2),size(actions_tpc,2));
            wlan(i).channel = a;   
            wlan(i).CCA = actions_cca(b);
            wlan(i).PTdBm = actions_tpc(c);  
        end
        tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
        % Compute the joint reward according to the average tpt obtained
        reward = mean(tpt_aux);
        progressive_tpt(r,iteration) = reward;
        progressive_fairness(r,iteration) = JainsFness(reward);
        % Assign reward to configurations set at each WLAN
        for i=1:N_WLANs
            reward_per_configuration(i,selected_arm(i)) = mean([reward_per_configuration(i,selected_arm(i)) reward]) + sqrt(2*log(iteration)/times_played_arm(i,selected_arm(i)));
            %regrets(i,iteration) = max([0 (max(tpt_aux)-reward)]);
        end
        avg_regret_evolution(iteration) = max(tpt_aux)-reward;%mean(regrets(:,iteration));
        % CHECK RESULTS AS THE NUMBER OF ITERATIONS INCREASES
%         if mod(iteration,T) == 0
%             % Apply the best configuration for each WLAN
%             for i=1:N_WLANs
%                 [~, index] = max(reward_per_configuration(i,:));
%                 [a,b,c] = val2indexes(index,NumChannels,size(actions_cca,2),size(actions_tpc,2));
%                 wlan(i).channel = a;
%                 wlan(i).CCA = actions_cca(b);
%                 wlan(i).PTdBm = actions_tpc(c);
%             end
%             % COMPUTE SNR AND TPT FOR EACH WLAN ACCORDING TO BACKOFF SIMULATION
%             tpt_aux = computeTpt(wlan,powMat,time_slots,noise); % bps
%             progressive_tpt(r,count) = mean(tpt_aux);
%             progressive_fairness(r,count) = JainsFness(tpt_aux);
%             count = count + 1;
%         end
         iteration = iteration + 1;
    end
    %%%%%%%%%%%%%%%%%%% MA-MAB implementation - END %%%%%%%%%%%%%%%%%%%
    
    % COMPUTE THE RESULTING THROUGHPUT AND FAIRNESS COEFFICIENT AT THIS ROUND
    % Apply the best configuration for each WLAN and compute the avg tpt and the fairness index
    for i=1:N_WLANs
        % Select configuration with minimum regret at each WLAN
        [~,index] = max(reward_per_configuration(i,:));
        [a,b,c] = val2indexes(states_wlan(index),NumChannels,size(actions_cca,2),size(actions_tpc,2));
        wlan(i).channel = a;
        wlan(i).CCA = actions_cca(b);
        wlan(i).PTdBm = actions_tpc(c);
    end
    tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
    avg_tpt_end(r) = mean(tpt_aux);
    % Check if results obtained are better than what we had at the beginning
    if avg_tpt_start(r) > 1.1*avg_tpt_end(r)
        %Go back to initial conf
        disp('Fail. Return to initial status')
        avg_tpt_end(r) = avg_tpt_start(r);
        fairness_end(r) = fairness_start(r);
    else
        fairness_end(r) = JainsFness(tpt_aux);
    end
    % Display new channels allocation
    for i=1:N_WLANs
        channels(i) = wlan(i).channel;       
    end
    for i=1:NumChannels
        sumChannels(i) = sum(channels==i);
    end
    disp('Channels selected per WLAN (at the end)')
    disp(channels)
    disp('Times a channel is occupied (at the end)')
    disp(sumChannels)
end

% Display average results
disp('MEAN IMPROVEMENTS IN OVERALL THROUGHPUT')
disp('-----------------------')
disp('AVG TPT START:')
disp(mean(avg_tpt_start))
disp('AVG TPT END:')
disp(mean(avg_tpt_end))
disp('MEAN IMPROVEMENTS IN FAIRNESS')
disp('-----------------------')
disp('AVG FAIRNESS START:')
disp(mean(fairness_start))
disp('AVG FAIRNESS END:')
disp(mean(fairness_end))

% Display evolution of the average throughput and fairness obtained
if R == 1
    figure
    progressive_tpt_smooth = smooth(progressive_tpt);
    plot(0:T:(size(progressive_tpt,2)*100)-1, progressive_tpt_smooth');
    xlabel('Number of iterations');
    ylabel('Average Throughput (bps)');
    title('Evolution of the throughput');

    figure
    progressive_fairness_smooth = smooth(progressive_fairness);
    plot(0:T:(size(progressive_tpt,2)*100)-1, progressive_fairness_smooth');
    hold on;
    xlabel('Number of iterations');
    ylabel('Jain''s Fairness Index');
    title('Evolution of the fairness');
else
    figure
    progressive_tpt_smooth = smooth(progressive_tpt);
    plot(0:T:(size(progressive_tpt,2)*100)-1, progressive_tpt_smooth);
    xlabel('Number of iterations');
    ylabel('Average Throughput (bps)');
    title('Evolution of the throughput');

    figure
    progressive_fairness_smooth = smooth(progressive_fairness);
    plot(0:T:(size(progressive_tpt,2)*100)-1, mean(progressive_fairness_smooth));
    hold on;
    xlabel('Number of iterations');
    ylabel('Jain''s Fairness Index');
    title('Evolution of the fairness');
end