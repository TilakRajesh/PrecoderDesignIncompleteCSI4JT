function SimulationMainLambdaGanesh(SNRdB_edge,threshold, BSs,  Nt,  UEs,  Nr)
datestr(now)
SNRdB_edge
threshold
rng('shuffle')
%% CVX setup
% cd cvx
% cvx_setup
% cd ..
%%
tic
MC = 1e3 % Number of iterations
% threshold = 3; % Relative Threshold [dB]
% save SNRdB_edge_Nt.mat SNRdB_edge BSs  Nt  UEs  Nr
SNRdB_edge = 15 % To avoid file read/write, change in this file CoMP_BS_Pwr.m
[Pmax,  ~, BSs, Nt, UEs, Nr, ~, ~] = CoMP_BS_Pwr();

C_ZF_FFB_FBH = zeros(MC, 1);
C_SSOCP_FFB_FBH = zeros(MC, 1);
C_SSOCP_LFB_LBH_PL_TriInEq = zeros(MC, 1);
C_SSOCP_LFB_LBH_PL = zeros(MC, 1);
C_SSOCP_LFB_LBH_PL_Dummy = zeros(MC, 1);
C_SSOCP_LFB_LBH = zeros(MC, 1);
C_PSO_FFB_FBH = zeros(MC, 1);
C_PSO_LFB_LBH = zeros(MC, 1);
PrecoderRate_SSOCP_FFB_FBH = zeros(MC, 1);
PrecoderRate_SSOCP_LFB_LBH_PL_TriInEq = zeros(MC,1);
PrecoderRate_SSOCP_LFB_LBH_PL = zeros(MC, 1);
PrecoderRate_SSOCP_LFB_LBH_PL_dummy = zeros(MC, 1);
PrecoderRate_SSOCP_LFB_LBH = zeros(MC, 1);
PrecoderRate_PSO_FFB_FBH = zeros(MC,1);
PrecoderRate_PSO_LFB_LBH = zeros(MC,1);


% C_SSOCP_FFB_FBH_tran = zeros(MC, 1);
% C_SSOCP_LFB_LBH_PL_tran = zeros(MC, 1);
% C_SSOCP_LFB_LBH_PL_Dummy_tran = zeros(MC, 1);
% C_SSOCP_LFB_LBH_tran = zeros(MC, 1);
% PrecoderRate_SSOCP_FFB_FBH_tran = zeros(MC, 1);
% PrecoderRate_SSOCP_LFB_LBH_PL_tran = zeros(MC, 1);
% PrecoderRate_SSOCP_LFB_LBH_PL_dummy_tran = zeros(MC, 1);
% PrecoderRate_SSOCP_LFB_LBH_tran = zeros(MC, 1);
%%
% matlabpool close force
% if (exist('local_scheduler_data','dir') == 0)
%     mkdir local_scheduler_data
% end
% sched=findResource('scheduler','type','local');
% sched.DataLocation = strcat(pwd,'/local_scheduler_data');
% matlabpool open local
% parfor k = 1:MC
for k = 1:MC
%     rng(1)
    [H,  PL] = CoMP_IID(BSs, Nt, UEs);
%     %% Make a diagonal matrix
%     [~,iMax]=max(PL,[],2);
%     Hnew = H;
%     PLnew = PL;
%     for iUe = 1:UEs
%         Hnew(iUe,iMax(iUe)) = H(iUe,iUe);
%         Hnew(iUe,iUe)= H(iUe,iMax(iUe));
%         PLnew(iUe,iMax(iUe)) = PL(iUe,iUe);
%         PLnew(iUe,iUe)= PL(iUe,iMax(iUe));
%     end
%     H = Hnew;
%     PL = PLnew;
    %%
    [activeSet] = RelativeThresholdBasedOnPathloss(BSs,  Nt,  UEs,  Nr,  PL,  threshold);
%     activeSet = eye(3);
    H_LFB = H .* activeSet; % Remove the inactive links
    %% ZF FFB+FBH
    W = H'*(H*H')^-1;
    W = NormalizeBF(W);
    W = LimitBSTransmitPower(W, Pmax);
    [C_ZF_FFB_FBH(k), ~, ~, ~, ~]  = CalculateRate(H, W, nan);
	%% OptionsFlag
    % OptionsFlag = 0; % No pathloss, use as is for Full and limited case
    % OptionsFlag = 1; % Use pathloss correctly (Antti's lambda)
    % OptionsFlag = 2; % Use pathloss in a dummy way (my stupid way)
%     %% Tran: WSRM: FFB+FBH
%     [W,PrecoderRate_SSOCP_FFB_FBH_tran(k)] = SolveFastConvAlgoWSRM(H, ones(size(activeSet)), PL, 0);
%     [C_SSOCP_FFB_FBH_tran(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %% Ganesh: SOCP: LFB + LBH with pathloss + TriInEq
    [W,PrecoderRate_SSOCP_LFB_LBH_PL_TriInEq(k)] = SolveFastConvAlgoWSRMLambdaGaneshLFB_LBH_PL_TriInEq(H_LFB, activeSet, PL, 1);
    [C_SSOCP_LFB_LBH_PL_TriInEq(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);        
%% Ganesh: SOCP: FFB + FBH: Tight upperbound
    [W,PrecoderRate_SSOCP_FFB_FBH(k)] = SolveFastConvAlgoWSRMLambdaGanesh(H, ones(size(activeSet)), PL, 0);
    [C_SSOCP_FFB_FBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
%     %% TRAN: SOCP: LFB + LBH
%     [W,PrecoderRate_SSOCP_LFB_LBH_tran(k)] = SolveFastConvAlgoWSRM(H_LFB, activeSet, PL, 0);
%     [C_SSOCP_LFB_LBH_tran(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %% Ganesh: SOCP: LFB + LBH
    [W,PrecoderRate_SSOCP_LFB_LBH(k)] = SolveFastConvAlgoWSRMLambdaGanesh(H_LFB, activeSet, PL, 0);
    [C_SSOCP_LFB_LBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
%     %% Tran: SOCP: LFB + LBH with pathloss correctly
%     [W,PrecoderRate_SSOCP_LFB_LBH_PL_tran(k)] = SolveFastConvAlgoWSRMLambda(H_LFB, activeSet, PL, 1);
%     [C_SSOCP_LFB_LBH_PL_tran(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);

    %% Ganesh: SOCP: LFB + LBH with pathloss
    [W,PrecoderRate_SSOCP_LFB_LBH_PL(k)] = SolveFastConvAlgoWSRMLambdaGanesh(H_LFB, activeSet, PL, 1);
    [C_SSOCP_LFB_LBH_PL(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
%     %% TRAN: SOCP: LFB + LBH with pathloss dummy approach
%     [W,PrecoderRate_SSOCP_LFB_LBH_PL_dummy_tran(k)] = SolveFastConvAlgoWSRM(H_LFB, activeSet, PL, 1);
%     [C_SSOCP_LFB_LBH_PL_Dummy_tran(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %% Ganesh: SOCP: LFB + LBH with pathloss dummy approach
    [W,PrecoderRate_SSOCP_LFB_LBH_PL_dummy(k)] = SolveFastConvAlgoWSRMLambdaGanesh(H_LFB, activeSet, PL, 2);
    [C_SSOCP_LFB_LBH_PL_Dummy(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
%     
    %% PSO: FFB+FBH
    [W,PrecoderRate_PSO_FFB_FBH(k)] = PSO_ParticleSwarmRetry(ones(size(activeSet)), H, Pmax);
    W = LimitBSTransmitPower(W, Pmax);
    [C_PSO_FFB_FBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %% PSO: LFB+LBH
    [W,PrecoderRate_PSO_LFB_LBH(k)] = PSO_ParticleSwarmRetry(activeSet, H_LFB, Pmax);
    W = LimitBSTransmitPower(W, Pmax);
    [C_PSO_LFB_LBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
end
toc
% eval(['save log_PSO_bugfix_Ganesh_celledgeSNR_',num2str(SNRdB_edge),'dB_thres_',num2str(threshold),'dB_NumBSAnt_',num2str(Nt),'.mat C_PSO_FFB_FBH C_PSO_LFB_LBH PrecoderRate_PSO_FFB_FBH PrecoderRate_PSO_LFB_LBH'])
eval(['save test_logAll_TriInEq_Ganesh_celledgeSNR_',num2str(SNRdB_edge),'dB_thres_',num2str(threshold),'dB_NumBSAnt_',num2str(Nt),'.mat C_ZF_FFB_FBH C_SSOCP_FFB_FBH C_SSOCP_LFB_LBH_PL C_SSOCP_LFB_LBH_PL_Dummy C_SSOCP_LFB_LBH C_PSO_FFB_FBH C_PSO_LFB_LBH PrecoderRate_SSOCP_FFB_FBH PrecoderRate_SSOCP_LFB_LBH_PL PrecoderRate_SSOCP_LFB_LBH_PL_dummy PrecoderRate_SSOCP_LFB_LBH PrecoderRate_PSO_FFB_FBH PrecoderRate_PSO_LFB_LBH  BSs Nt UEs Nr PrecoderRate_SSOCP_LFB_LBH_PL_TriInEq C_SSOCP_LFB_LBH_PL_TriInEq'])
end
