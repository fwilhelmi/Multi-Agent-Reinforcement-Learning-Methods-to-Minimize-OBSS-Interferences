% MASTER THESIS - Minimizing OBSS Interferences using Q-learning and Correlated Equilibria
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

% I: Show coexistence effects on throughput 

clc
clear all

N = 1000; % Number of iterations for computing the average
N_WLANs = 10; % Maximum number of WLANs
%NumChannels = [3 6 9]; % Maximum number of channels available
NumChannels = [1 6 11];
TxPower = -20:10:20; % Tx power (dBm)
CCA = -30:10:-10; % Clear Channel Assesment value (dBm)
B = 20e6; %Bandwidth per channel (bps)
noise = -100; %Noise (dBm)
AverageTpt = zeros(size(NumChannels,2),N_WLANs); %Average throughput for each pair (NumChannels,N_WLANs)
AggregateTpt = zeros(size(NumChannels,2),N_WLANs); %Aggregate throughput for each pair (NumChannels,N_WLANs
colors = [1 1 0;1 0 0;0 1 0;0 0 1;0 1 1; 1 0 1];

%CSMA/CA information
CW = 10;
BO=floor(rand()*2^CW);

approach = 2;

time_slots = 100;

total_size = [size(NumChannels,2) size(TxPower,2) size(CCA,2)];

%% FIRST APPROACH
if approach == 1
    % REPEAT CALCULATIONS FOR EACH NUMBER OF CHANNELS
    for c=1:size(NumChannels,2)
        avgTpt = zeros(N_WLANs,N); %
        aggTpt = zeros(N_WLANs,N);
        % PERFORM CALCUATIONS "N" TIMES AND TAKE THE AVERAGE
        for n=1:N
            % REPEAT CALCULATIONS FOR EACH NUMBER OF WLANs
            for k=1:N_WLANs
                wlan = GenerateNetwork3D(k, NumChannels(c), B);
                powMat = PowerMatrix(wlan); % dBm
                tpt_aux = computeTpt(wlan,powMat,time_slots,noise); 
                aggTpt(k,n) = aggTpt(k,n) + mean(tpt_aux);
                avgTpt(k,n) = aggTpt(k,n)/k;                   
            end
        end
        for k=1:N_WLANs
            mean_avg(c,k)=mean(avgTpt(k,:));
            std_avg(c,k)=std(avgTpt(k,:));
            mean_agg(c,k)=mean(aggTpt(k,:));
            std_agg(c,k)=std(aggTpt(k,:));
            AggregateTpt(c,k) = (sum(aggTpt(k,:))/N); 
            AverageTpt(c,k) = (sum(avgTpt(k,:))/N); 
        end
    end
    l = {}; %legend
    for i=1:size(NumChannels,2)
        l = [l; ['NumChannels = ' num2str(NumChannels(i))]];
    end       
%     disp('Average Tpt/WLAN')
%     disp(AverageTpt)

%% SECOND APPROACH: modifying TPC
elseif approach == 2
    NumChannels = 6;
    % REPEAT CALCULATIONS FOR EACH CCA VALUE
    for c=1:size(TxPower,2)
        avgTpt = zeros(N_WLANs,N); %
        aggTpt = zeros(N_WLANs,N);
        % PERFORM CALCUATIONS "N" TIMES AND TAKE THE AVERAGE
        for n=1:N
            % REPEAT CALCULATIONS FOR EACH NUMBER OF WLANs
            for k=1:N_WLANs
                wlan = GenerateNetwork3D(k, NumChannels, B);
                for x=1:k, wlan(x).PTdBm = TxPower(c); end                    
                powMat = PowerMatrix(wlan); % dBm
                tpt_aux = computeTpt(wlan,powMat,time_slots,noise); 
                aggTpt(k,n) = aggTpt(k,n) + mean(tpt_aux);
                avgTpt(k,n) = aggTpt(k,n)/k;                   
            end
        end
        for k=1:N_WLANs
            mean_avg(c,k)=mean(avgTpt(k,:));
            std_avg(c,k)=std(avgTpt(k,:));
            mean_agg(c,k)=mean(aggTpt(k,:));
            std_agg(c,k)=std(aggTpt(k,:));
            AggregateTpt(c,k) = (sum(aggTpt(k,:))/N); 
            AverageTpt(c,k) = (sum(avgTpt(k,:))/N); 
        end
    end
    l = {}; %legend
    for i=1:size(TxPower,2)
        l = [l; ['TxPower = ' num2str(TxPower(i))]];
    end
    disp('Average Tpt/WLAN')
    disp(AverageTpt)
    
%% THIRD APPROACH: modifying CCA
elseif approach == 3
    NumChannels = 6;
    % REPEAT CALCULATIONS FOR EACH CCA VALUE
    for c=1:size(CCA,2)
        avgTpt = zeros(N_WLANs,N); %
        aggTpt = zeros(N_WLANs,N);
        % PERFORM CALCUATIONS "N" TIMES AND TAKE THE AVERAGE
        for n=1:N
            % REPEAT CALCULATIONS FOR EACH NUMBER OF WLANs
            for k=1:N_WLANs                
                wlan = GenerateNetwork3D(k, NumChannels, B);
                for x=1:k, wlan(x).CCA = CCA(c); end                 
                powMat = PowerMatrix(wlan); % dBm
                tpt_aux = computeTpt(wlan,powMat,time_slots,noise); 
                aggTpt(k,n) = aggTpt(k,n) + mean(tpt_aux);
                avgTpt(k,n) = aggTpt(k,n)/k;                   
            end
        end
        for k=1:N_WLANs
            mean_avg(c,k)=mean(avgTpt(k,:));
            std_avg(c,k)=std(avgTpt(k,:));
            mean_agg(c,k)=mean(aggTpt(k,:));
            std_agg(c,k)=std(aggTpt(k,:));
            AggregateTpt(c,k) = (sum(aggTpt(k,:))/N); 
            AverageTpt(c,k) = (sum(avgTpt(k,:))/N); 
        end
    end
    l = {}; %legend
    for i=1:size(CCA,2)
        l = [l; ['CCA = ' num2str(CCA(i))]];
    end  
    disp('Average Tpt/WLAN')
    disp(AverageTpt)
end

%% PLOT THE RESULTS
for i=1:total_size(approach)
        subplot(2,2,1)
        plot(1:N_WLANs, smooth(AverageTpt(i,:)));  
        hold on;        
        subplot(2,2,2)
        plot(1:N_WLANs, smooth(AggregateTpt(i,:)));
        hold on;
        subplot(2,2,3)
        errorbar(1:N_WLANs,smooth(mean_avg(i,:)),std_avg(i,:),'ko-','LineWidth',1,'Markersize',5,'Color',colors(i,:)); 
        hold on;
        subplot(2,2,4)
        errorbar(1:N_WLANs,smooth(mean_agg(i,:)),std_agg(i,:),'ko-','LineWidth',1,'Markersize',5,'Color',colors(i,:)); 
        hold on;
end 
subplot(2,2,1)
xlabel('NWLANs');
ylabel('Average Capacity (bps)');
title('Average Capacity vs #WLANs');
legend(l);
axis([1 N_WLANs 0 max(max(AverageTpt))]);

subplot(2,2,2)
xlabel('NWLANs');
ylabel('Aggregate Capacity (bps)');
title('Aggregate Capacity vs #WLANs')
legend(l);
axis([1 N_WLANs 0 max(max(AggregateTpt))]);

subplot(2,2,3)
xlabel('NWLANs');
ylabel('Average Capacity (bps)');
title('Average Capacity vs #WLANs with error');
legend(l);
axis([1 N_WLANs 0 max(max(AverageTpt))]);

subplot(2,2,4)
xlabel('NWLANs');
ylabel('Aggregate Capacity (bps)');
title('Aggregate Capacity vs #WLANs with error')
legend(l);
axis([1 N_WLANs 0 max(max(AggregateTpt))]);