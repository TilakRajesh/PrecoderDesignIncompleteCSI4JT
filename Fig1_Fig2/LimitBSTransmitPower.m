function [WOut] = LimitBSTransmitPower(W,Pmax)

% [Pmax, ~, BSs, Nt, ~, ~, sumPwrFlag, perAntennaPwrFlag] = CoMP_BS_Pwr();
[~, ~, BSs, Nt, ~, ~, sumPwrFlag, perAntennaPwrFlag] = CoMP_BS_Pwr();

if perAntennaPwrFlag == 1
    maxPwrTemp = zeros(1,BSs*Nt);
    for k = 1:BSs*Nt
        maxPwrTemp(k) = norm(W(k,:),'fro')^2;
    end
else
    maxPwrTemp = zeros(1,BSs);
    for k = 1:BSs
        maxPwrTemp(k) = norm(W((k-1)*Nt+1:k*Nt,:),'fro')^2;
    end
end
if sumPwrFlag
    maxPwr = max(maxPwrTemp);
    factor = Pmax / maxPwr;
else
    %Power per-BS
    factor = Pmax ./ maxPwrTemp; 
    % Take care of the divide by zero, if a BS is not serving any UE
    factor(isinf(factor)) = 0; 
    factor  = repmat(factor',1,3); 
end
 
WOut = sqrt(factor) .* W;
