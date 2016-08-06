function [WoutRet, maxRateAchievedSOCP] = SolveFastConvAlgoWSRMLambdaGanesh(H, activeSet, PL, OptionsFlag)
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
%%
Wout = zeros(BSs*Nt, UEs);
WoutRet = Wout;
maxRateAchievedSOCP = 0;
W_AS = activeSet';
numBFCoeff = sum(sum(activeSet)); % CVX uses this
maxTries = 5;

for randBF = 1:maxTries % Make 5 attempts
    %% Random BF initialization
    W_init = sqrt(0.5)*(randn(size(H_rand_init)) + 1j*randn(size(H_rand_init)));
    W_init = W_init.*activeSet';
    W_init = LimitBSTransmitPower(W_init, Pmax);
    [~, SINR_init, ~, ~, hw] = CalculateRate(H_rand_init, W_init, N0);
    pk0 = real(hw); qk0 = imag(hw); % Method-2: Q-WSRM
    t_k_n = 1 + SINR_init;
    betak0 = ((pk0.^2 + qk0.^2)./(t_k_n - 1));
    %%
    n=1;
    alpha_k = ones(UEs,1);
    maxIter = 200;
    objectiveValueWSR = zeros(1,maxIter);
    SINR=zeros(UEs, maxIter);
    WSRM= zeros(1,maxIter);
    while (n <= maxIter && (n<=6 ||  (objectiveValueWSR(n-1)-objectiveValueWSR(n-6))/objectiveValueWSR(n-1)>eps))
        %% Solve
        cvx_begin
        cvx_solver sdpt3%Gurobi%mosek% %
        cvx_precision best
        cvx_quiet(true);
        variable W(numBFCoeff, 1) complex; % vec of [BSs*Nt x UEs]
        variable t_k(UEs,1); % [UEs x 1]
        variable beta_k(UEs,1); % [UEs x 1]
        expression eachNormTerm(UEs,1)
        expression pk(UEs,1);
        expression qk(UEs,1);
        
        maximize(geo_mean(t_k))
        subject to
        
        %%
        ind = 1;
        for ue=1:UEs
            allSigUeOnly = 0;
            for bs=1:BSs*Nt
                if W_AS(bs,ue) == 1
                    allSigUeOnly = allSigUeOnly + H(ue,bs)*W(ind);
                    ind = ind + 1;
                end
            end
            ind2 = 1;  intUEVector = []; flag = 0; flag_PL_lambda = 0;
            for intUEs=1:UEs
                allIntUes = 0;
                allIntUes_PL_lambda = [];
                for bs=1:BSs*Nt
                    if W_AS(bs,intUEs) == 1
                        if intUEs ~= ue
                            switch OptionsFlag
                                case 0
                                    allIntUes = allIntUes + H(ue,bs)*W(ind2);
                                    flag = 1;
                                case 1
                                    if H(ue,bs) == 0
                                        % Use the pathloss value correctly
                                        allIntUes_PL_lambda = [allIntUes_PL_lambda; PL_Nt(ue,bs)*W(ind2)];
                                        flag_PL_lambda = 1;
                                    else
                                        allIntUes = allIntUes + H(ue,bs)*W(ind2);
                                        flag = 1;
                                    end
                                case 2
                                    % In this dummy approach, the PL values
                                    % are already added to H matrix
                                    allIntUes = allIntUes + H(ue,bs)*W(ind2);
                                    flag = 1;
                            end
                        end
                        ind2 = ind2 + 1;
                    end
                end
                if (flag == 1)
                    intUEVector = [intUEVector; allIntUes];
                    flag = 0;
                end
                if (flag_PL_lambda == 1)
                    intUEVector = [intUEVector; allIntUes_PL_lambda];
                    flag_PL_lambda = 0;
                end
            end
            intVector= [sqrt(N0);intUEVector; ];
            eachNormTerm(ue) = allSigUeOnly;
            pk(ue) = real(allSigUeOnly);
            qk(ue) = imag(allSigUeOnly);
            intVector= [2 * intVector; (beta_k(ue)-1)];
            {intVector,(beta_k(ue)+1)} == complex_lorentz(numel(intVector));
            %             (t_k_n(ue)^(1/alpha_k(ue))) + (1/alpha_k(ue))*(t_k_n(ue)^((1/alpha_k(ue))-1))*(t_k(ue) - t_k_n(ue)) <= ...
            %                 + (2*(pk0(ue)/betak0(ue))*(pk(ue)-pk0(ue))) ...
            %                 + (2*(qk0(ue)/betak0(ue))*(qk(ue)-qk0(ue)))...
            %                 + ((pk0(ue)^2+qk0(ue)^2)/betak0(ue))*(1-((beta_k(ue)-betak0(ue))/betak0(ue))) ...
            %                 + 1 ;
            t_k(ue) <= ...
                + (2*(pk0(ue)/betak0(ue))*(pk(ue)-pk0(ue))) ...
                + (2*(qk0(ue)/betak0(ue))*(qk(ue)-qk0(ue)))...
                + ((pk0(ue)^2+qk0(ue)^2)/betak0(ue))*(1-((beta_k(ue)-betak0(ue))/betak0(ue))) ...
                + 1 ;
        end
        % Per-BS Power constraint
        for i=1:BSs*Nt
            WperBS = [];
            for perBs=i:BSs*Nt:UEs*BSs*Nt;
                if W_AS(perBs) == 1
                    index = sum(W_AS(1:perBs));
                    WperBS = [WperBS; W(index)];
                end
            end
            if ~isempty(WperBS)
                {WperBS, sqrt(Pmax)} == complex_lorentz(numel(WperBS));
            end
        end
        cvx_end
        %% Update
        if (strfind(cvx_status,'Solved'))
            t_k_n = t_k;
            betak0 = beta_k;
            pk0 = pk;
            qk0 = qk;
            SINR(:,n)=abs(eachNormTerm).^2./beta_k.^2;
            objectiveValueWSR(n)=sum(log2(1+SINR(:,n)));
            WSRM(n) = sum(log2(t_k));
            if WSRM(n) > maxRateAchievedSOCP
                % Store the maximum value
                maxRateAchievedSOCP = WSRM(n);
                %% Demapp the W
                k=1;
                for ue=1:UEs
                    for perBs=1:BSs*Nt
                        if activeSet(ue,perBs) == 1
                            Wout(perBs,ue) = W(k);
                            k=k+1;
                        end
                    end
                end
                WoutRet = Wout;
            end
        else
            cvx_status
            randBF
            H
            n
            c = num2str(fix(clock));
            c(ismember(c,' ,.:;!')) = [];
            eval(['save Failed_SolveFastConvAlgoWSRMLambdaGanesh',num2str(c),'.mat'])
            break;
        end
        n = n+1;
    end
end
end
