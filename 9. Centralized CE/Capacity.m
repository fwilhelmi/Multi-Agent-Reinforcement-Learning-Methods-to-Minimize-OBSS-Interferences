function C = Capacity(B, sinr)
%%Capacity - Returns the theoretical capacity given a bandwidth and a SNR
    % B - Available Bandwidth (Hz) 
    % sinr - Signal to Interference plus Noise Ratio (-)
    C = B * log2(1+sinr);
end

% Theoretical Capacity
% B = [5e6 10e6 15e6 20e6];
% for j = 1:size(B,2)
%     for i=0:100
%         C(j,i+1) = Capacity(B(j), 10^(i/10));
%     end
% end
% for i=1:size(B,2)
%     plot(0:100,C(i,:))
%     hold on
% end
% legend('5MHz','10MHz','15MHz','20MHz')
% xlabel('E/N0 (dB)')
% ylabel('Capacity (bps)')
% title('Theoretical Shannon Capacity')