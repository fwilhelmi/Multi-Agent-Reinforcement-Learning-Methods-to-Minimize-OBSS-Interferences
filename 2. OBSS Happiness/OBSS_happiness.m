% MASTER THESIS - Multi-Agent Reinforcement Learning Methods to Minimize OBSS Interferences
% Author - Francesc Wilhelmi
% Tutors - Boris Bellalta & Anders Jonsson

% II: Show that enhancements are possible and how hard are to be achieved

clc
clear all

approach = 2; %approach to be selected (1 or 2)
plot = 0;

N = 10; %Number of repetitions
N_WLANs = 1:10; %Maximum number of WLANs
NumChannels = [1 6 11]; %Number of channels
B = 20e6; %Bandwidth required (Hz)
CCA = -82; %CCA level by default (dBm)
noise = -100; %Noise leve dBm
maxIter = 50; %Maximum convergence time

%CSMA/CA information
time_slots = 100;

%% FIRST APPROACH: 
% WLANs change the channel randomly if not "happy"
convergenceTime = zeros(size(NumChannels,2), size(N_WLANs,2));

stop = false;

for c=1:size(NumChannels,2)
    maxIt = zeros(size(N_WLANs,2), N);       
    for n=1:N
        disp('-----------')
        txt = ['ROUND ' num2str(n) '/' num2str(N)];
        disp(txt)
        for k=2:size(N_WLANs,2)
            wlan = GenerateNetwork3D(N_WLANs(k),NumChannels(c),B);          
            found = false;
            it = 1;
            force_stop = false;
            while it < maxIter + 1
                %total_tpt = sum(tpt_aux);
                %minimum_condition = max(tpt_aux)*0.2/(total_tpt);
                for i=1:N_WLANs(k)
                    powMat = PowerMatrix(wlan); % dBm
                    tpt_aux = computeTpt(wlan,powMat,time_slots,noise);
                    minimum_condition = 0.3*Capacity(wlan(i).BW,db2pow(150))*min([1,NumChannels(c)/N_WLANs(k)]);
                    % EVALUATE IF WLAN IS "HAPPY"
                    if tpt_aux(i) < minimum_condition % Check if happiness condition is accomplished by each WLAN                        
                        found = true;
                        %if rand() < 1/N_WLANs(k)
                        if approach == 1
                            if NumChannels(c) == 1
                                force_stop = true;
                                break; 
                            end
                            change = true;
                            while change
                                newCh = ceil(NumChannels(c)*rand());
                                if (wlan(i).channel ~= newCh)
                                    wlan(i).channel = newCh;
                                    change = false;
                                end
                            end
                        elseif approach == 2
                            for i=1:N_WLANs
                                channels(i) = wlan(i).channel;       
                            end
                            for i=1:NumChannels
                                channels_occupancy(i) = sum(channels==i);
                            end
                            [~,index] = min(channels_occupancy);
                            wlan(i).channel = index;
                            
%                             powerChannel = PowerPerChannel(wlan,NumChannels(c),powMat);
%                             [val, ~] = min(powerChannel(i,:));
%                             if ~powerChannel(i,wlan(i).channel)==val
%                                 indexes = (powerChannel(i,:)==val);
%                                 indexes(wlan(i).channel) = 0;
%                                 feasible_options = [];
%                                 for l=1:size(indexes,2)
%                                     if indexes(l)~=0, feasible_options = [feasible_options l]; end
%                                 end
%                                 [~,newCh] = datasample(feasible_options,1);
%                                 if pow2db(powerChannel(i,wlan(i).channel)) > pow2db(powerChannel(i,newCh)), wlan(i).channel = newCh; end
%                             end
                        else
                            disp('Approach introduced does not exist!');
                        end
                    end
                end                
                maxIt(k,n) = it;
                % If not found, we are done
                if ~found, break;
                else it = it + 1; end  
            end
            if NumChannels(c) == 1 && force_stop, maxIt(k,n) = maxIter; end
        end
    end
    % Compute the average results obtained during simulations
    for k=1:size(N_WLANs,2)
        convergenceTime(c,k) = sum(maxIt(k,:))/N; 
        %if convergenceTime(c,k) == maxIter || (k>1 && convergenceTime(c,k-1)==1e10), convergenceTime(c,k) = Inf; end
    end
end
% Plot results
if plot
    N_WLANs = 1:10;
    NumChannels = [1 6 11];
    maxIter = 50;
    figure
    l = {}; %legend
    for i=1:size(NumChannels,2)
        plot(N_WLANs, convergenceTime(i,:),'-*');
        hold on;
        %plot(N_WLANs, convTimeSmart(i,:),'--x');
        l = [l; ['Random approach (' num2str(NumChannels(i)) ' channels)']];
        %l = [l; ['Smart approach (' num2str(NumChannels(i)) ' channels)']];
    end
    xlabel('NWLANs');
    ylabel('Convergence Time (iterations)');
    title('Reaching a compliance state on an OBSS');
    legend(l);
    axis([1 max(N_WLANs) 0 maxIter]);
end