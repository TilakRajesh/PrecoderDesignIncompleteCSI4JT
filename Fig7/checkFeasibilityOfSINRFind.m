function [feasibleFlag,Wout] = checkFeasibilityOfSINRFind(H,gammaSINRUE,activeSet,N0,whichFuncCalled,OptionsFlag,PL_Nt)

[Pmax, ~, BSs, Nt, UEs, ~, ~, ~] = CoMP_BS_Pwr();
Wout = zeros(BSs*Nt,UEs);
numBFCoeff = sum(sum(activeSet));
W_AS = activeSet';
cvx_begin
cvx_quiet( true )
cvx_solver gurobi %sdpt3
variable W(numBFCoeff, 1) complex; % vec of [BSs*Nt x UEs]
expression eachNormTerm(UEs,1);
find W;
subject to
%%
ind = 1;
for ue=1:UEs
    if gammaSINRUE(ue) ~= 0
        allSigUeOnly = 0;
        for bs=1:BSs*Nt
            if W_AS(bs,ue) == 1
                allSigUeOnly = allSigUeOnly + H(ue,bs)*W(ind);
                ind = ind + 1;
            end
        end
        ind2 = 1;  intUEVector = []; flag = 0; allVector = []; flag_PL_lambda = 0;
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
%                 intUEVector = [intUEVector; allIntUes_PL_lambda];
                intUEVector = [intUEVector; sqrt(numel(allIntUes_PL_lambda))*allIntUes_PL_lambda];
                flag_PL_lambda = 0;
            end
        end
        eachNormTerm(ue) = allSigUeOnly;
        allVector= [eachNormTerm(ue); intUEVector; sqrt(N0)]; % Note that the find algorithm requires all that is needed for a given UE
        % SOC contraints
        {allVector, sqrt(1+(1/gammaSINRUE(ue)))*eachNormTerm(ue)} == complex_lorentz(numel(allVector));
        %% hard code the constraints for testing
        %             %Get all the things related to one UE
        %             allVecPerUE = [];
        %             switch ue
        %                 case 1
        %                     W_temp = W(1:3); W_tempInt2 = W(4:6); W_tempInt3 = W(7:9);
        %                     allVecPerUE = [H(ue,:)*W_temp; H(ue,:)*W_tempInt2; H(ue,:)*W_tempInt3];
        %                 case 2
        %                     W_temp = W(4:6); W_tempInt1 = W(1:3); W_tempInt3 = W(7:9);
        %                     allVecPerUE = [H(ue,:)*W_tempInt1; H(ue,:)*W_temp; H(ue,:)*W_tempInt3];
        %                 case 3
        %                     W_temp = W(7:9); W_tempInt1 = W(1:3); W_tempInt2 = W(4:6);
        %                     allVecPerUE = [H(ue,:)*W_tempInt1; H(ue,:)*W_tempInt2; H(ue,:)*W_temp];
        %             end
        %             allVecPerUE = [allVecPerUE; sqrt(N0)];
        %             {allVecPerUE, sqrt(1+(1/gammaSINRUE(ue)))*H(ue,:)*W_temp} == complex_lorentz(numel(allVecPerUE));
    end
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

if ~isempty(strfind(cvx_status,'Solved'))
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
    feasibleFlag = 1;
else
    % Infeasible
    feasibleFlag = 0;
    if whichFuncCalled
        disp('Fail: Trouble')
        c = num2str(fix(clock));
        c(ismember(c,' ,.:;!')) = [];
        % Its okay to be infeasible, so we don't go to that hyperrectangle
%         eval(['save BB_caller_CVX_checkWhy',num2str(c),'.mat'])
    end
end

end
