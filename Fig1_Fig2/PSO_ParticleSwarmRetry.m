function [WoutRet,maxRateAchievedPSO] = PSO_ParticleSwarmRetry(activeSet, H, Pmax)
maxRetrires = 5;
globalWBest = zeros(size(H'));
maxRateAchievedPSO = 0;
for iMaxRetryIter=1:maxRetrires
    %% Configure the necessary zeros
    W = activeSet';
    %% Settings
    numberOfIterations = 500;
    numberOfParticles = 30; % 20 to 40
    numberOfVariables = sum(sum(W))*2; % 2 is due to complex numbers, here the real part and the imaginary parts are optimized individually
    xmax = 1/max(max(abs(H)));
    xmin = -xmax;
    alpha = 1;
    deltaT = 1;
    % Minimization, hence choose a large value
    evaluatedBestMinimumOfEachParticle = Inf*ones(numberOfParticles,1);
    evaluatedGlobalBestMinimum = Inf;
    bestPositionOfEachParticle = zeros(numberOfParticles,numberOfVariables);
    globalBestPosition = 0;
    c1 = 2;
    c2 = 2;
    maxInertiaWeight = 1.1;
    minInertiaWeight = 0.4;
    inertiaWeight = maxInertiaWeight;
    beta = 0.99;
    vmax = (xmax - xmin)/deltaT;
    bestMinEveryPSOIteration = nan*ones(1,numberOfIterations);
    %% Initialize Positions and Velocities
    [position, velocity] = PSO_InitializeParticles(numberOfParticles, numberOfVariables, xmax, xmin, alpha, deltaT);
    %% RPSO
    % sIdx = 0; % Based on Algorithm 15.7 in the book "Fundamentals of Computational swarm Intelligence
    % tau = 0;
    % tr = 30;
    %% Evaluate Particle and update the best position and global position
    for iteration = 1:numberOfIterations
        %% RPSO
        %     if (sIdx ~= tau)
        %         [position(sIdx,:), velocity(sIdx,:)] = PSO_InitializeParticles(1, numberOfVariables, xmax, xmin, alpha, deltaT);
        %     end
        %     sIdx = mod(sIdx,tr)+1;
        %     if iteration > 100
        %         % decay tr
        %         tr = tr-1;
        %         if (tr < 10)
        %             tr = 10;
        %         end
        %     end
        %%
        for i = 1:numberOfParticles
            [evaluatedMinimum, WBest] = PSO_EvaluateParticle(H, W, position(i,:), Pmax);
            if evaluatedMinimum < evaluatedBestMinimumOfEachParticle(i)
                evaluatedBestMinimumOfEachParticle(i) = evaluatedMinimum;
                bestPositionOfEachParticle(i,:) = position(i,:);
            end
            if evaluatedMinimum < evaluatedGlobalBestMinimum
                evaluatedGlobalBestMinimum = evaluatedMinimum;
                bestMinEveryPSOIteration(iteration) = evaluatedGlobalBestMinimum;
                globalBestPosition = position(i,:);
                globalWBest = WBest;
            else
                bestMinEveryPSOIteration(iteration) = evaluatedGlobalBestMinimum;
            end
        end
        %% Update particle velocities and positions
        
        for i = 1:numberOfParticles
            for j = 1:numberOfVariables
                %             if numel(globalBestPosition) == 1
                %                 globalBestPosition
                %             end
                velocity(i,j) = inertiaWeight*velocity(i,j) + c1*rand*(bestPositionOfEachParticle(i,j) - position(i,j))/deltaT ...
                    + c2*rand*(globalBestPosition(j) - position(i,j))/deltaT;
                if abs(velocity(i,j)) > vmax
                    if velocity(i,j) > 0
                        velocity(i,j) = vmax; % restrict velocity
                    else
                        velocity(i,j) = -vmax; % restrict velocity
                    end
                end
                position(i,j) = position(i,j) + velocity(i,j) * deltaT;
            end
        end
        
        if inertiaWeight < minInertiaWeight
            inertiaWeight = minInertiaWeight;
        else
            inertiaWeight = inertiaWeight*beta;
        end
        
    end
    %evaluatedBestMinimumOfEachParticle, evaluatedGlobalBestMinimum%, globalWBest
    %%
    globalWBest = LimitBSTransmitPower(globalWBest, Pmax);
    [rateAchievedPSOTemp, ~, ~, ~,~] = CalculateRate(H, globalWBest, nan);
    if rateAchievedPSOTemp > maxRateAchievedPSO
        % Store the maximum value
        maxRateAchievedPSO = rateAchievedPSOTemp;
        WoutRet = globalWBest;
    end
end
end
