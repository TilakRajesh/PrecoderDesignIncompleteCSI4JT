function RetW = NormalizeBF(W)
% NumUE = size(W,2);
% RetW = zeros(size(W));
% for ue=1:NumUE
%     RetW(:,ue) = W(:,ue)/norm(W(:,ue));
% end
%%
% if size(W,1) == 1
%     size(W)
% end
NumUE = size(W,2);
RetW = zeros(size(W));
for ue=1:NumUE
    RetW(:,ue) = W(:,ue)/norm(W(:,ue));
end

end