% MASTER THESIS - A ML Solution For OBSS Coexistence in Dense Environments
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

clc
clear all

% Global variables for the problem
N_WLANs = 5;
NumChannels = 6;
noise = -100;
B = 20e6;

N = 1000; %Number of iterations
epsilon = 0.1;
gamma = 0.8;
alpha = 0.7;
% Definition of actions:
%   - 11 possible channels
%   - 2 different levels of CCA
%   - 2 different levels of TxPower
actions_ch = 1:NumChannels;
actions_cca = [-82 -72]; %dBm
actions_tpc = [15 20]; %dBm
%actions_bw = [20e6 40e6];

% Each state represents an [i,j,k] combination for indexes on "channels", "CCA" and "TxPower"
states_wlan = 1:(size(actions_ch,2)*size(actions_cca,2)*size(actions_tpc,2));
% The total states matrix is given by the configurations at each WLAN
states_total = allcomb(states_wlan,states_wlan,states_wlan,states_wlan,states_wlan);
reward_per_state = zeros(1,size(states_total,1));

R=1;

avg_tpt_start = zeros(1,R);
fairness_start = zeros(1,R);
avg_tpt_end = zeros(1,R);
fairness_end = zeros(1,R);

time_slots = 100;

for r=1:R
    disp('-----------')
    txt = ['ROUND ' num2str(r) '/' num2str(R)];
    disp(txt)
    time_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    disp(time_start)
    disp('-----------')
    % Generate network of N WLANs
    %wlan = GenerateNetwork3D(N_WLANs, NumChannels, B);
    load('wlan.mat');
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
    
    powMat = PowerMatrix(wlan); %                  
     
    tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
    avg_tpt_start(r)=mean(tpt_aux);
    fairness_start(r)=JainsFness(tpt_aux);

    for k=1:size(states_total,1)  
        selected_action = states_total(k,:);
        for i=1:N_WLANs
            % Choose action according to Q-val            
            % Change parameters according to the action obtained
            [a b c]=val2indexes(selected_action(i),size(actions_ch,2),size(actions_cca,2),size(actions_tpc,2));
            wlan(i).channel = a;
            wlan(i).CCA = actions_cca(b);
            wlan(i).PTdBm = actions_tpc(c);
        end
%         val_curr_state=find(states_total(:,1)==curr_state(1) & states_total(:,2)==curr_state(2) & states_total(:,3)==curr_state(3) & states_total(:,4)==curr_state(4) & states_total(:,5)==curr_state(5));
%         val_next_state=find(states_total(:,1)==next_state(1) & states_total(:,2)==next_state(2) & states_total(:,3)==next_state(3) & states_total(:,4)==next_state(4) & states_total(:,5)==next_state(5));
        powMat = PowerMatrix(wlan); %
        tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
        reward_per_state(k) = mean(tpt_aux)*JainsFness(tpt_aux);        
    end
    %%%%%%%%%%%%%%%%%%% NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [val,indx]=max(reward_per_state);
    final_states=states_total(indx,:);
    for i=1:N_WLANs
        [a b c]=val2indexes(final_states(i),size(actions_ch,2),size(actions_cca,2),size(actions_tpc,2));
        wlan(i).channel=a;
        wlan(i).CCA=actions_cca(b);
        wlan(i).PTdBm=actions_tpc(c);
    end
    powMat = PowerMatrix(wlan); %
    tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
    avg_tpt_end(r)=mean(tpt_aux);
    fairness_end(r)=JainsFness(tpt_aux);
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