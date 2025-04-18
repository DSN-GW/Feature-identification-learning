%% sweepmodel
% Complex model, one free parameter. Taken from Rosenblau et al., 2018
% ER = y * (ER + a*PE) + (1-y) * OP

%% Function
function [sweepVec, sumSquareError] = sweepmodel(data, param, sweepCorr,startVal)
    learnVec = zeros(length(data),1);
    sweepVec = zeros(length(data),1);
    resetVal = data(:,3);
    %global startVal
    
    for i = 1:length(data)
        % For each new profile you need to reset
        if(resetVal(i) == 1)
        if (length(param) ~= 2)
        spaceVec = zeros( length(sweepCorr), 1) +  startVal;
        elseif(length(param) == 2); spaceVec = zeros( length(sweepCorr), 1) +  param(2); end
        end
        
        learnVec(i) = spaceVec(data(i,10));           
        delta = data(i,5) - learnVec(i);

        spaceVec = spaceVec + param(1) * delta * sweepCorr(:, data(i,10));

        sweepVec(i) = learnVec(i);
    end
    % Add bounding that is lacking in fminsearch
    if param(1) < 0 || param(1) > 1
        sumSquareError = 10^25;
    else
        sumSquareError = sum( ( data(:,4) - sweepVec ).^2 );
    end
end