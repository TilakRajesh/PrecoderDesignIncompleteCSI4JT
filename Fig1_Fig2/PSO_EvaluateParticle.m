function [error, WBest] = PSO_EvaluateParticle(H, WIn, position, Pmax)
users = size(H,1);
I = eye(users);
bsNt = size(H,2);
WBest = zeros(size(WIn));
%% Form the W matrix where ones are present based on the 2 consecutive positions
n = 1;
for k = 1:bsNt
    for m = 1:users
        if (WIn(k,m) == 1)
            WBest(k,m) = position(n) + 1i * position(n+1);
            n = n + 2;
        end
    end
end
%% 
WBest = LimitBSTransmitPower(WBest, Pmax);
[SumRate, ~, ~,~,~] = CalculateRate(H, WBest, nan);
%Object function
error = -SumRate;
% error = -min(SINR);
%%
end
