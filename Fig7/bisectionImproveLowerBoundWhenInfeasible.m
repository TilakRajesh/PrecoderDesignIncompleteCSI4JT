function [gammaStar] = bisectionImproveLowerBoundWhenInfeasible(H,Q,activeSet,N0,OptionsFlag,PL_Nt)
% gammaStar is the feasible gamma in G, so
[~,~, ~,~, UEs, ~, ~, ~] = CoMP_BS_Pwr();
e = eye(UEs); gammaStar = zeros(UEs,1);
gammaMin = Q(:,2);
epsilon = 1e-2;
% parfor iUE=1:UEs
for iUE=1:UEs
    a =  gammaMin + (Q(iUE,1) - Q(iUE,2))*e(:,iUE);
    
    L = gammaMin; U = a;
    [feasibilityFlag, ~] = checkFeasibilityOfSINRFind(H,a,activeSet,N0,0,OptionsFlag,PL_Nt);
    if feasibilityFlag
        gammaStar(iUE) = a(iUE);
        continue;
    end
    while (norm(U-L,2) > epsilon)
        T = (L+U)/2;
        [feasibilityFlag, ~] = checkFeasibilityOfSINRFind(H,T,activeSet,N0,0,OptionsFlag,PL_Nt);
        if feasibilityFlag == 1
            L = T;
        else
            U = T;
        end
    end
   gammaStar(iUE) = U(iUE);
end

end