% MASTER THESIS - A ML Solution For OBSS Coexistence in Dense Environments
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

clc
clear all

N_WLANs = 5;
NumChannels = 6;
noise = -100; % dBm
B = 20e6;

R = 1;

possible_ch = 1:NumChannels;
possible_cca = [-82 -70]; %dBm
possible_tpc = [15 25]; %dBm
% Set of possible configurations (strategies that each WLAN can implement)
Configurations = zeros(size(possible_ch,2)*size(possible_cca,2)*size(possible_tpc,2),3);
index = 1;
a = 1;
while(a<size(possible_ch,2)+1)
    b = 1;
    while(b<size(possible_cca,2)+1)
        c = 1;
        while(c<size(possible_tpc,2)+1)
            Configurations(index,1) = possible_ch(a);
            Configurations(index,2) = possible_cca(b);
            Configurations(index,3) = possible_tpc(c);
            index = index + 1;
            c = c + 1;
        end
        b = b + 1;
    end
    a = a + 1;
end

% Create a data structure with all the possible configurations for each WLAN
c={};
for i=1:N_WLANs
    c{i} = 1:size(Configurations,1);
end
AllCombinations = allcomb(c{:});
% Remove "bad" combinations (i.e. all WLANs choose the same channel)
% ...
% ...
time_slots = 100;

avg_tpt_start = zeros(1,R);
avg_tpt_end = zeros(1,R);
fairness_start = zeros(1,R);
fairness_end = zeros(1,R);

for r=1:R
    disp('-----------')
    txt = ['ROUND ' num2str(r) '/' num2str(R)];
    disp(txt)
    disp('-----------')
    % Generate map of WLANs
    %wlan = GenerateNetwork3D(N_WLANs, NumChannels, B);
    load('wlan.mat');
    powMat = PowerMatrix(wlan);
    tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
    % Resulting throughput
    avg_tpt_start(r) = mean(tpt_aux);
    disp('Average Tpt (Start)')
    disp(avg_tpt_start(r))
    fairness_start(r) = JainsFness(tpt_aux);
    disp('Fairness (Start)')    
    disp(fairness_start(r))

    %DrawNetwork3D(wlan,powMat);
    % Possible actions to be done
    
    % Compute payoffs and best response
    % The payoff matrix contains the average throughput if the Jain's 
    % coefficient is higher than 0.8. Otherwise, it is NaN. 
    % 
    %   Example of payoff Matrix (2 WLANs and 3 actions)
    %       - An action is a set of configurations for the Channel, the CCA
    %       and the TPC

    %                --------------------------------------
    %               | Ch = 1;    | Ch = 2;    | Ch = 3;    |
    %   WLAN1/WLAN2 | CCA = -82; | CCA = -82; | CCA = -82; |
    %               | TPC = 20;  | TPC = 20;  | TPC = 20;  |
    %   ---------------------------------------------------
    %  | Ch = 1;    |            |            |            |
    %  | CCA = -82; |   3e8 bps  |   5e8 bps  |   7e8 bps  |
    %  | TPC = 20;  |            |            |            |
    %   ---------------------------------------------------
    %  | Ch = 2;    |            |            |            |           
    %  | CCA = -82; |   5e8 bps  |   3e8 bps  |   5e8 bps  | 
    %  | TPC = 20;  |            |            |            |
    %   ---------------------------------------------------
    %  | Ch = 3;    |            |            |            |           
    %  | CCA = -82; |   7e8 bps  |   5e8 bps  |   3e8 bps  | 
    %  | TPC = 20;  |            |            |            |
    %   ---------------------------------------------------
    disp('Computing payoff matrix...');
    payoff = zeros(1,size(Configurations,1)^N_WLANs);    
    aux = 1;
    for i=1:size(AllCombinations)
        % Keep track on the computation status
        if mod(i,ceil(size(AllCombinations,1)/10))==0
            disp([num2str(10*aux) '% done']);
            aux = aux + 1;
        end
        if length(unique(Configurations(AllCombinations(i,:),1))) < length(Configurations(AllCombinations(i,:),1))
            payoff(i)=NaN;
            continue;
        end
        for n=1:N_WLANs
            wlan(n).channel = Configurations(AllCombinations(i,n),1);
            wlan(n).CCA = Configurations(AllCombinations(i,n),2);
            wlan(n).PTdBm = Configurations(AllCombinations(i,n),3);
        end
        powMat = PowerMatrix(wlan);
        tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
        if JainsFness(tpt_aux)>0.8, payoff(i) = mean(tpt_aux);
        else payoff(i) = NaN; end    
    end

    disp('Finding the best configuration...');
    [val,index] = max(payoff);
    found = AllCombinations(index,:);
    bestConf={};
    for i=1:N_WLANs
        bestConf{i} = Configurations(found(i),:);
    end
    % Assign to each WLAN the computed configuration
    disp('Best Configuration found:')
    for i=1:N_WLANs
        disp('--------')
        txt = ['WLAN #' num2str(i)];
        disp(txt)
        wlan(i).channel = bestConf{i}(1);    
        wlan(i).CCA = bestConf{i}(2);
        wlan(i).PTdBm = bestConf{i}(3);
        txt2 = ['Channel = '  num2str(bestConf{i}(1)) '; CCA = '  num2str(bestConf{i}(2)) '; TPC = '  num2str(bestConf{i}(3))];
        disp(txt2)
    end

    % Resulting throughput
    powMat = PowerMatrix(wlan);
    tpt_end = computeTpt(wlan,powMat,time_slots,noise);
    avg_tpt_end(r) = mean(tpt_end);
    disp('Average Tpt (End)')
    disp(avg_tpt_end(r))
    fairness_end(r) = JainsFness(tpt_end);
    disp('Fairness (End)')    
    disp(fairness_end(r))
    % Draw the network after configuration is done
    %DrawNetwork3D(wlan,powMat);

    %Distance between WLANs
    d = zeros(N_WLANs,N_WLANs);
    for i=1:N_WLANs
        for j=1:N_WLANs
            if i~=j
                d(i,j) = sqrt(sum(([wlan(i).x wlan(i).y wlan(i).z] - [wlan(j).x wlan(j).y wlan(j).z]).^2));
            end
        end
    end
end