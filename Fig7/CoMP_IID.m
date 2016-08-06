function [H,PL]=CoMP_IID(BSs,Nt,UEs)
R = 500; % Cell Radius
Dij = CalculateDistancesBSUE(BSs, UEs, R);
H = zeros(UEs,BSs*Nt);
PL = zeros(UEs,BSs);
for user=1:UEs
    for bs = 1:BSs
        % Shadow fading
        fading=8*randn;  % dB
        fading=10^(-fading/10);
        % TR 36.942 Release 10
        pathloss = 128.1 + 37.6*log10(Dij(user,bs)*10^(-3)); % R is in kms
        pathloss=10.^(-pathloss/10);      % lineal

        G = (1/sqrt(2)) * (randn(1,Nt) + 1i * randn(1,Nt));

%         G_tx=9; % dB
%         G_tx=10^(G_tx/10);  % lineal
        G_tx=1;

        % Correlation
        rho=0.5;
        C=rho*ones(Nt,Nt);
        C=C-diag(diag(C))+eye(Nt);

        H(user,(bs-1)*Nt+1:bs*Nt)=G*C^(1/2)*sqrt(G_tx*pathloss*fading);
        PL(user,bs) = sqrt(G_tx*pathloss*fading); % RSSI known to the transmitter
    end
end
end
