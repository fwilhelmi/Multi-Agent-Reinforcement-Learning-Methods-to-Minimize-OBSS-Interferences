% MASTER THESIS - Multi-Agent Reinforcement Learning Methods to Minimize OBSS Interferences
% Decentralized Q-learning - SHARED Approach
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

clc
clear all

disp('DECENTRALIZED Q-LEARNING - SHARED APPROACH')
disp('-----------------------')
% Global variables for the problem
N_WLANs = 5;
NumChannels = 6; %Number of available channels (from 1 to NumChannels)
noise = -100; % Floor noise (dBm)
B = 20e6; % Available bandwidth (Hz)

N = 1000; % Number of iterations for Q-learning
epsilon = 0.1; % Exploration coefficient
gamma = 0.95; % Discount rate
alpha = 0.3; % Learning rate
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
%states_total = allcomb(states_wlan,states_wlan,states_wlan,states_wlan,states_wlan);

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
progressive_tpt = zeros(R,N/T);
progressive_fairness = zeros(R,N/T);

for r=1:R % We repeat the experiment during 'R' repetitions for taking the average 
    disp('-----------')
    txt = ['ROUND ' num2str(r) '/' num2str(R)];
    disp(txt)
    time_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    disp(time_start)
    disp('-----------')
    % Generate network of N WLANs and measure Tpt and Fairness
    disp('Creating map of OBSS...')
    %wlan = GenerateNetwork3D(N_WLANs, NumChannels, B);
    load('wlan.mat');
    %load(data_wlan)
    % wlan = data_wlan.wlan;
    powMat = PowerMatrix(wlan); %
    
    count = 1;
    
    %%%%%%%%%%%%%%%%%%% Q-learning implementation - START %%%%%%%%%%%%%%%%%
    
    % Initialize initial and next state arrays of each WLAN 
    % (each state represents a configuration [channel,CCA,TPC])
    initial_state = zeros(1,N_WLANs); 
    current_state = zeros(1,N_WLANs); 
    next_state = zeros(1,N_WLANs);
    % Define the initial state of each WLAN
    for i=1:N_WLANs    
        [~,index_cca] = find(actions_cca==wlan(i).CCA);
        [~,index_tpc] = find(actions_tpc==wlan(i).PTdBm);
        initial_state_indexes = [wlan(i).channel index_cca index_tpc]; %Take indexes of each selected value    
        initial_state(i) = indexes2val(initial_state_indexes(1),initial_state_indexes(2),initial_state_indexes(3),size(actions_ch,2),size(actions_cca,2)); %Convert indexes into a numeric state
    end
    curr_state=initial_state; 
    % Compute initial tpt and fairness
    tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
    avg_tpt_start(r) = mean(tpt_aux);
    fairness_start(r) = JainsFness(tpt_aux);
    % Initiliaze Q-Values Table Q(s), which contains the value for performing an action on each state
    % At step 0 we can have several states (e.g. WLAN #1 has selected channel 1, transmits at 20 dBm and its CCA is -82)
    Qval = {};
    for i=1:N_WLANs
        Qval{i}=zeros(1,size(states_wlan,2));
    end
    % Apply the algorithm in parallel for each WLAN
    throughput = zeros(N_WLANs, N); % Keep track of the tpt on each WLAN for each iteration
    iteration = 1;  
    powMat = PowerMatrix(wlan); %
    reward = zeros(1,N);
    while(iteration <= N)
        for i=1:N_WLANs            
            % IMPLEMENT Q-LEARNING AT EACH WLAN
            % Choose action according to the current Q-val of each WLAN
            selected_action = ChooseAction(Qval{i},curr_state(i),actions_ch,actions_cca,actions_tpc,epsilon);
            % Change parameters according to the action obtained
            wlan(i).channel = selected_action(1);   
            wlan(i).CCA = actions_cca(selected_action(2));
            wlan(i).PTdBm = actions_tpc(selected_action(3));
            % Prepare the next state according to the actions performed on the current state
            [~,index_cca] = find(actions_cca==wlan(i).CCA);
            [~,index_tpc] = find(actions_tpc==wlan(i).PTdBm);
            next_state(i) = indexes2val(wlan(i).channel,index_cca,index_tpc,size(actions_ch,2),size(actions_cca,2));
        end
        % Compute the reward with the average throughput obtained (consider fairness)
        throughput(:,iteration) = computeTpt(wlan,powMat,time_slots,noise); % bps
        reward(iteration) = mean(throughput(:,iteration))*JainsFness(throughput(:,iteration)');
        for i=1:N_WLANs
            Qval{i}(curr_state(i))=(1-alpha)*Qval{i}(curr_state(i)) + alpha*(reward(iteration)+gamma*Qval{i}(next_state(i))); %Update Q
        end
        curr_state = next_state;
        % CHECK RESULTS AS THE NUMBER OF ITERATIONS INCREASES
        if mod(iteration,T) == 0
            % Apply the best configuration for each WLAN
            for i=1:N_WLANs
                [~, index] = max(Qval{i});
                [a,b,c] = val2indexes(index,NumChannels,size(actions_cca,2),size(actions_tpc,2));
                wlan(i).channel = a;
                wlan(i).CCA = actions_cca(b);
                wlan(i).PTdBm = actions_tpc(c);
            end
            % COMPUTE SNR AND TPT FOR EACH WLAN ACCORDING TO BACKOFF SIMULATION
            tpt_aux = computeTpt(wlan,powMat,time_slots,noise); % bps
            progressive_tpt(r,count) = mean(tpt_aux);
            progressive_fairness(r,count) = JainsFness(tpt_aux);
            count = count + 1;
        end
        iteration = iteration + 1;
    end
    
    %%%%%%%%%%%%%%%%%%% Q-learning implementation - END %%%%%%%%%%%%%%%%%%%
    
    % COMPUTE THE RESULTING THROUGHPUT AND FAIRNESS COEFFICIENT
    % Apply the best configuration for each WLAN
    for i=1:N_WLANs
        [~, index] = max(Qval{i});
        [a,b,c] = val2indexes(index,NumChannels,size(actions_cca,2),size(actions_tpc,2));
        wlan(i).channel = a;
        wlan(i).CCA = actions_cca(b);
        wlan(i).PTdBm = actions_tpc(c);
    end
    % COMPUTE SNR AND TPT FOR EACH WLAN ACCORDING TO BACKOFF SIMULATION
    tpt_aux = computeTpt(wlan,powMat,time_slots,noise); % bps
    % Check if results obtained are better than what we had at the beginning    
    avg_tpt_end(r) = mean(tpt_aux);
    if avg_tpt_start(r) > 1.1*avg_tpt_end(r)
        %Go back to initial conf
        disp('Fail. Return to initial status')
        avg_tpt_end(r) = avg_tpt_start(r);
        fairness_end(r) = fairness_start(r);
    else
        % Keep the computed value and calculate fairness
        fairness_end(r) = JainsFness(tpt_aux);
    end
    % Display new channels allocation
    disp('Channels selected per WLAN (at the end)')
    for i=1:N_WLANs
        channels(i) = wlan(i).channel;       
    end
    for i=1:NumChannels
        sumChannels(i) = sum(channels==i);
    end
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
    %progressive_tpt_smooth = smooth(progressive_tpt);
    plot(0:T:(size(progressive_tpt,2)*100)-1, mean(progressive_tpt));
    xlabel('Number of iterations');
    ylabel('Average Throughput (bps)');
    title('Evolution of the throughput');

    figure
    %progressive_fairness_smooth = smooth(progressive_fairness);
    plot(0:T:(size(progressive_tpt,2)*100)-1, mean(progressive_fairness));
    hold on;
    xlabel('Number of iterations');
    ylabel('Jain''s Fairness Index');
    title('Evolution of the fairness');
end