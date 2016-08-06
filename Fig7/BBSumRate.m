function [WOptimal,U,L, Usave, Lsave] = BBSumRate(H,activeSet,PL,maxIterBB,OptionsFlag)
[Pmax, N0, BSs, Nt, UEs, ~, ~, ~] = CoMP_BS_Pwr();
PLfactor = min(min(PL));
%% OptionsFlag
% OptionsFlag = 0; % No pathloss, use as is for Full and limited case
% OptionsFlag = 1; % Use pathloss correctly (Antti's lambda)
% OptionsFlag = 2; % Use pathloss in a dummy way (my stupid way)
switch OptionsFlag
    case 0
        % Do nothing
        H_rand_init = H;
        PL_Nt = nan; % just to make sure that other parts of the code work =)
    case 1
        % Use pathloss correctly
        PL_Nt = zeros(UEs, BSs*Nt);
        for bs=1:BSs
            % Just repeat
            PL_Nt(:,(bs-1)*Nt+1:bs*Nt) = repmat(PL(:,bs),1,Nt);
        end
        PL_Nt = PL_Nt.*~activeSet; % Take the pathloss values where they are zeros
        H_rand_init = H + PL_Nt.*~activeSet; % used for random BF init.
        PL_Nt = PL_Nt/PLfactor;
    case 2
        % Dummy approach
        PL_Nt = zeros(UEs, BSs*Nt);
        for bs=1:BSs
            % Just repeat
            PL_Nt(:,(bs-1)*Nt+1:bs*Nt) = repmat(PL(:,bs),1,Nt);
        end
        H = H + PL_Nt.*~activeSet; % used for random BF init.
        H_rand_init = H;
        PL_Nt = PL_Nt/PLfactor;
    otherwise
        disp('Unknown OptionsFlag')
        return;
end
%% overwrite to help convergence problems when N0=4e-15
H = H/PLfactor;
H_rand_init = H_rand_init/PLfactor;
N0=N0/PLfactor^2;
%% Branch and Bound from Algorithm 1 in paper WSRMax for MISO DL in cellular networks via branch and bound
epsilon = 1e-1;
% maxIterBB = 1e3%100; 
Usave = zeros(maxIterBB,1); Lsave = zeros(maxIterBB,1);
% utopia = abs(diag(H*H').^2)/N0; 
% wMRT = functionMRT(H);
% utopia = abs(diag(H*wMRT*sqrt(Pmax)).^2)/N0; % Seems correct but poor results
% If there is only one active user then let all the BSs transmit to this one user, therefore, there is no interference
utopia = abs(diag(H*H'))*BSs*Pmax/N0;

% utopia = 1e5*ones(UEs,1); %test
Qcurr = [utopia 1e-5*ones(UEs,1)];% [UpperVec LowerVec]; All rectangles are in this format
[feasibilityFlag, ~] = checkFeasibilityOfSINRFind(H,Qcurr(:,2),activeSet,N0,1,OptionsFlag,PL_Nt);
[Qcurr(:,1)] = bisectionImproveLowerBoundWhenInfeasible(H,Qcurr,activeSet,N0,OptionsFlag,PL_Nt);
if feasibilityFlag == 1
    U = -EvaluateBound(Qcurr(:,2));
    L = -EvaluateBound(Qcurr(:,1));
else
    U = 0;
    L = 0;
    disp('Fail: Trouble')
    c = num2str(fix(clock));
    c(ismember(c,' ,.:;!')) = [];
    eval(['save FeasibilityFailedInBegin',num2str(c),'.mat'])
end
Bk{1}.Q = Qcurr;
Bk{1}.U = U;
Bk{1}.L = L;
iter = 1;
Usave(iter) = U; Lsave(iter) = L; 
while ((U - L > epsilon) && (iter <= maxIterBB)) %|| ((~feasibilityFlagQ1) && (~feasibilityFlagQ2)) %|| (sum((Qcurr(:,1)-Qcurr(:,2) > epsilonSINR*ones(UEs,1))))
%     U
%     L
    iter = iter + 1;
    %Branching 
    indSizeBk = size(Bk,2);
    saveIndex = 0;
    for indBk = 1:indSizeBk
        if L == Bk{indBk}.L;
            Qcurr = Bk{indBk}.Q;
            saveIndex = indBk;
            break;
        end
    end

    %Split the longest edgeHar
    [~,index] = max(Qcurr(:,1)-Qcurr(:,2));
    mid = (Qcurr(index,1) - Qcurr(index,2))/2;
    Q1 = Qcurr;
    Q1(index,2) = Qcurr(index,2) + mid;
    Q2 = Qcurr;
    Q2(index,1) = Qcurr(index,2) + mid;
    
    % Perform Bisection for the lower bound
    [Q1(:,1)] = bisectionImproveLowerBoundWhenInfeasible(H,Q1,activeSet,N0,OptionsFlag,PL_Nt);
    %Bounding for Q1
    [feasibilityFlagQ1, ~] = checkFeasibilityOfSINRFind(H,Q1(:,2),activeSet,N0,1,OptionsFlag,PL_Nt);
    if feasibilityFlagQ1 == 1
        U_Q1 = -EvaluateBound(Q1(:,2));
        L_Q1 = -EvaluateBound(Q1(:,1));
    else
        U_Q1 = 0; 
        L_Q1 = 0;
    end
    Bk{saveIndex}.Q = Q1; Bk{saveIndex}.L = L_Q1; Bk{saveIndex}.U = U_Q1; % step 3.c [16]
    %Bounding for Q2
    [Q2(:,1)] = bisectionImproveLowerBoundWhenInfeasible(H,Q2,activeSet,N0,OptionsFlag,PL_Nt);
    [feasibilityFlagQ2, ~] = checkFeasibilityOfSINRFind(H,Q2(:,2),activeSet,N0,1,OptionsFlag,PL_Nt);
    if ~feasibilityFlagQ1 || ~feasibilityFlagQ2
        disp('Fail: How can Q1 or Q2 be infeasible? something kaput!')
        c = num2str(fix(clock));
        c(ismember(c,' ,.:;!')) = [];
        % Its okay, as the hyperrectangle will be going back to the
        % previous bigger one. %Confirm this analysis.
%         eval(['save Q1_or_Q2_is_infeasible_checkWhy_',num2str(c),'.mat'])
    end
    if feasibilityFlagQ2 == 1
        U_Q2 = -EvaluateBound(Q2(:,2));
        L_Q2 = -EvaluateBound(Q2(:,1));
    else
        U_Q2 = 0;
        L_Q2 = 0;
    end
    Bk{indSizeBk+1}.Q = Q2; Bk{indSizeBk+1}.U = U_Q2; Bk{indSizeBk+1}.L = L_Q2; % step 3.c [16]
    for indBk = 1:size(Bk,2)
        UkAll(indBk) = Bk{indBk}.U;
        LkAll(indBk) = Bk{indBk}.L;
    end
    U = min(UkAll);
    L = min(LkAll);
    Usave(iter) = U; Lsave(iter) = L; % Saving logs
end

%% Calculate the beamformer
if feasibilityFlagQ1
%     sumRateBB_upper_Q1 = sum(log2(1+Q1(:,2)))
%     sumRateBB_lower_Q1 = sum(log2(1+Q1(:,1)))
    [~, WOptimal] = checkFeasibilityOfSINRFind(H,Q1(:,2),activeSet,N0,1,OptionsFlag,PL_Nt);
else
        disp('Fail: How can gammaMin of Q1 be infeasible?')
end
if feasibilityFlagQ2
%     sumRateBB_upper_Q2 = sum(log2(1+Q2(:,2)))
%     sumRateBB_lower_Q2 = sum(log2(1+Q2(:,1)))
    [~, WOptimal] = checkFeasibilityOfSINRFind(H,Q2(:,2),activeSet,N0,1,OptionsFlag,PL_Nt);
else
    disp('Fail: How can gammaMin of Q2 be infeasible?')
end
% sumRateBB_upper = max([sumRateBB_upper_Q1 sumRateBB_upper_Q2]);
% sumRateBB_lower = max([sumRateBB_lower_Q1 sumRateBB_lower_Q2]);

% [rate_BB, sinr_BB, ~, ~, ~] = CalculateRate(H,WOptimal, N0);
% figure
% plot(Usave(1:iter-1),'r--');hold on
% plot(Lsave(1:iter-1),'b--');  xlabel('Bound iterations'); ylabel(' min. (-SumRate)')
% [W,maxRateAchievedSOCP] = SolveFastConvAlgoWSRM(H, activeSet, PL, 0 );
% [rate_BB, sinr_BB, ~, ~, ~] = CalculateRate(H,WOptimal, nan);
% legend('Upper Bound','Lower Bound','','location','best')
%%
% Lsave
% save forPlot_1Iteration.mat Usave Lsave 