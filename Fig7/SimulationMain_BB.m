function SimulationMain_BB(SNRdB_edge,threshold)
datestr(now)
rng('shuffle')
% %% CVX setup
% cd cvx
% cvx_setup
% cd ..
%%
tic
convergence = 0;
if convergence
    MC = 4% enable=1 for getting the convergence plots
    maxIterBB = 1e3;
else
    MC = 1e2
    maxIterBB = 1e3;
end
[Pmax, ~, BSs, Nt, UEs, Nr, ~, ~] = CoMP_BS_Pwr();
maxRateAchievedSOCP_FFB_FBH = zeros(MC, 1);
C_SSOCP_FFB_FBH = zeros(MC, 1);
C_BB_FFB_FBH = zeros(MC, 1);
C_BB_FFB_FBH_UB = zeros(MC, 1);
C_BB_FFB_FBH_LB = zeros(MC, 1);

maxRateAchievedSOCP_LFB_LBH = zeros(MC, 1);
C_SSOCP_LFB_LBH = zeros(MC, 1);
C_BB_LFB_LBH = zeros(MC, 1);
C_BB_LFB_LBH_UB = zeros(MC, 1);
C_BB_LFB_LBH_LB = zeros(MC, 1);

maxRateAchievedSOCP_LFB_LBH_PL = zeros(MC, 1);
C_SSOCP_LFB_LBH_PL = zeros(MC, 1);
C_BB_LFB_LBH_PL_UB = zeros(MC, 1);
C_BB_LFB_LBH_PL_LB = zeros(MC, 1);
C_BB_LFB_LBH_PL = zeros(MC, 1);

C_ZF_FFB_FBH = zeros(MC, 1);
% C_PSO_FFB_FBH = zeros(MC, 1);
%matlabpool close 
c = num2str(fix(clock));
c(ismember(c,' ,.:;!')) = [];
% foldername = ['local_scheduler_data',num2str(c)];
% mkdir(foldername)
% pause(5)
% % sched=findResource('scheduler','type','local');
% sched.DataLocation = fullfile(pwd,foldername);
% matlabpool open local
parfor k = 1:MC
% for k = 1:MC
    rng('shuffle')
    [H,  PL] = CoMP_IID(BSs, Nt, UEs);
    [activeSet] = RelativeThresholdBasedOnPathloss(BSs,  Nt,  UEs,  Nr,  PL,  threshold);
    H_LFB = H .*activeSet;
    %% ZF FFB+FBH
    W = H'*(H*H')^-1;
    W = NormalizeBF(W);
    W = LimitBSTransmitPower(W, Pmax);
    [C_ZF_FFB_FBH(k), ~, ~, ~, ~]  = CalculateRate(H, W, nan);
    %% OptionsFlag
    % OptionsFlag = 0; % No pathloss, use as is for Full and limited case
    % OptionsFlag = 1; % Use pathloss correctly (Antti's lambda)
    % OptionsFlag = 2; % Use pathloss in a dummy way (my stupid way)
    %% SOCP: FFB + FBH: Tight upperbound
    OptionsFlag = 0;
    [W,maxRateAchievedSOCP_FFB_FBH(k)] = SolveFastConvAlgoWSRMLambdaGanesh(H, ones(size(activeSet)), PL, OptionsFlag);
    [C_SSOCP_FFB_FBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %  BB: FFB + FBH: Tight the ultimate upperbound
    if convergence
        [W,C_BB_FFB_FBH_UB(k),C_BB_FFB_FBH_LB(k), Usave_FFB_FBH(k,:), Lsave_FFB_FBH(k,:)] = BBSumRate(H,ones(size(activeSet)),PL,maxIterBB,OptionsFlag);
    else
        [W,C_BB_FFB_FBH_UB(k),C_BB_FFB_FBH_LB(k), ~, ~] = BBSumRate(H,ones(size(activeSet)),PL,maxIterBB,OptionsFlag);
    end
    [C_BB_FFB_FBH(k), ~, ~, ~, ~]= CalculateRate(H, W, nan);
    %% SOCP: LFB + LBH: Tight upperbound
    OptionsFlag = 0;
    [W,maxRateAchievedSOCP_LFB_LBH(k)] = SolveFastConvAlgoWSRMLambdaGanesh(H_LFB, activeSet, PL, OptionsFlag);
    [C_SSOCP_LFB_LBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %  BB: LFB + LBH: Tight the ultimate upperbound
    if convergence
        [W,C_BB_LFB_LBH_UB(k),C_BB_LFB_LBH_LB(k), Usave_LFB_LBH(k,:), Lsave_LFB_LBH(k,:)] = BBSumRate(H_LFB,activeSet,PL,maxIterBB,OptionsFlag);
    else
        [W,C_BB_LFB_LBH_UB(k),C_BB_LFB_LBH_LB(k), ~, ~] = BBSumRate(H_LFB,activeSet,PL,maxIterBB,OptionsFlag);
    end
    [C_BB_LFB_LBH(k), ~, ~, ~, ~]= CalculateRate(H, W, nan);
    %% SOCP: LFB + LBH: Pathloss including correctly
    OptionsFlag = 1;
    [W,maxRateAchievedSOCP_LFB_LBH_PL(k)] = SolveFastConvAlgoWSRMLambdaGaneshLFB_LBH_PL_TriInEq(H_LFB, activeSet, PL, OptionsFlag);
    [C_SSOCP_LFB_LBH_PL(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
    %  BB: LFB + LBH: Tight the ultimate upperbound
    if convergence
        [W,C_BB_LFB_LBH_PL_UB(k),C_BB_LFB_LBH_PL_LB(k), Usave_LFB_LBH_PL(k,:), Lsave_LFB_LBH_PL(k,:)] = BBSumRate(H_LFB,activeSet,PL,maxIterBB,OptionsFlag);
    else
        [W,C_BB_LFB_LBH_PL_UB(k),C_BB_LFB_LBH_PL_LB(k), ~,~] = BBSumRate(H_LFB,activeSet,PL,maxIterBB,OptionsFlag);
    end
    [C_BB_LFB_LBH_PL(k), ~, ~, ~, ~]= CalculateRate(H, W, nan);
    %% For plotting convergence
%     if convergence
%         save forPlotFFB_FBH_LFB_LBH_PL_1Iteration_longRune1e3.mat ...
%             Usave_FFB_FBH Lsave_FFB_FBH maxRateAchievedSOCP_FFB_FBH ...
%             Usave_LFB_LBH Lsave_LFB_LBH maxRateAchievedSOCP_LFB_LBH ...
%             Usave_LFB_LBH_PL Lsave_LFB_LBH_PL maxRateAchievedSOCP_LFB_LBH_PL
%         break;
%     end
    
%     %% PSO: FFB+FBH
%     [~, W] = PSO_ParticleSwarm(ones(size(activeSet)), H, Pmax);
%     W = LimitBSTransmitPower(W, Pmax);
%     [C_PSO_FFB_FBH(k), ~, ~, ~, ~] = CalculateRate(H, W, nan);
end
toc
eval(['save log_sdpt3_Ganesh_BB_MC_',num2str(MC),'_time_',num2str(c),'_celledgeSNR_',num2str(SNRdB_edge),'dB_thres_',num2str(threshold),'dB_NumBSAnt_',num2str(Nt),'.mat C_SSOCP_FFB_FBH C_BB_FFB_FBH C_BB_FFB_FBH_UB C_BB_FFB_FBH_LB C_SSOCP_LFB_LBH C_BB_LFB_LBH C_BB_LFB_LBH_UB C_BB_LFB_LBH_LB     C_SSOCP_LFB_LBH_PL C_BB_LFB_LBH_PL_UB C_BB_LFB_LBH_PL_LB C_BB_LFB_LBH_PL C_ZF_FFB_FBH BSs Nt UEs Nr Pmax convergence maxRateAchievedSOCP_FFB_FBH maxRateAchievedSOCP_LFB_LBH maxRateAchievedSOCP_LFB_LBH_PL'])

end
