function [activaHNum] = RelativeThresholdBasedOnPathloss(BSs,  Nt,  UEs,  Nr,  PL,  threshold)

PL_Nt = zeros(UEs, BSs*Nt);
for bs=1:BSs
    PL_Nt(:,(bs-1)*Nt+1:bs*Nt) = repmat(PL(:,bs),1,Nt);
end
HdB = 10*log10(PL_Nt.^2);

%%
% threshold = 10;%[0 3 5 10 15 20 40];
len = size(threshold,2);
activaHNum = zeros(UEs,BSs*Nt,len); % encapsulate all the activaH1-6
%%
for u=1:len
    activaH=ones(UEs,BSs*Nt);
    for i=1:UEs
        maxH=max(HdB(i,:));
        % Compares the current maximum value and attenuation of each channel MISO
        for k=1:BSs*Nt
            if (maxH-HdB(i,k))>threshold(u)
                % In this case there is a spending power to reach mobile
                % Do not need. Therefore the channel OFF
                activaH(i,k)=0;
            end
        end
    end
    activaHNum(:,:,u)=activaH; % Save all in one matrix :)
end
activaHNum = squeeze(activaHNum);
end
