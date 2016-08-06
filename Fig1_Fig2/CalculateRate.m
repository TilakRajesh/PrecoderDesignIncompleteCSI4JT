function [rate, sinr, sig, intf, diagB] = CalculateRate(H,W, N0)
if isnan(N0)
    N0 = CoMP_Rx_NoiseInit();
end
B = H*W;
diagB = diag(B);
NdiagB = B - diag(diag(B));
sig = abs(diagB).^2;
intf = (sum((abs(NdiagB)).^2,2)+N0);
sinr = sig./intf;
% if isnan(sinr)
%     sinr
% end
rate=sum(log2(1+sinr));

end
