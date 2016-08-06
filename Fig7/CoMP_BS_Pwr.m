function [Pdmax, N0, BSs, Nt, UEs, Nr, sumPwrFlag, perAntennaPwrFlag] = CoMP_BS_Pwr()

N0 = CoMP_Rx_NoiseInit();
R=500; % in meters
% SNRdB_edge = 15; % dB
% load SNRdB_edge_Nt.mat SNRdB_edge BSs Nt UEs Nr % To speed up: comment this
SNRdB_edge = 15; BSs = 3; Nt = 1; UEs = 3; Nr = 1;
% TR 36.942 Release 10
pathloss = 128.1 + 37.6*log10(R*10^(-3)); % R is in kms
PdB = SNRdB_edge+(pathloss+10*log10(N0));
Pdmax = 10.^(PdB/10); % in Watts
% 1= % This only means at least one of the BSs is transmitting at maximum power
% 0= Per BS power constraint % NEVER used yet!
sumPwrFlag = 1; % This only means at least one of the BSs is transmitting at maximum power
%% Per-Antenna power constriant, 
% 1= Enable
% 0= Disable
perAntennaPwrFlag = 1;
end
