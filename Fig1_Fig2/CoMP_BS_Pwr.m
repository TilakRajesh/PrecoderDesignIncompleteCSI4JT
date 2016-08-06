function [Pdmax, N0, BSs, Nt, UEs, Nr, sumPwrFlag, perAntennaPwrFlag] = CoMP_BS_Pwr()

N0 = CoMP_Rx_NoiseInit();
R=500; % in meters
%%
SNRdB_edge = 15; % dB
BSs=3; Nt=1; UEs=3; Nr=1;
% load SNRdB_edge_Nt.mat SNRdB_edge BSs Nt UEs Nr
%%
% TR 36.942 Release 10
pathloss = 128.1 + 37.6*log10(R*10^(-3)); % R is in kms
PdB = SNRdB_edge+(pathloss+10*log10(N0));
Pdmax = 10.^(PdB/10); % in Watts
% 1= Sum power constriant 
% 0= Per BS power constraint
sumPwrFlag = 1; 
%% Per-Antenna power constriant, 
% 1= Enable
% 0= Disable
perAntennaPwrFlag = 1;
end
