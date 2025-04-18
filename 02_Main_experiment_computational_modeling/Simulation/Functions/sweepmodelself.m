%% fitting_fminsearch_sweepmodelself
% Complex model, one free parameter. Taken from Rosenblau et al., 2018
% ER = y * (ER + a*PE) + (1-y) * OP
%
% Output is a vector with [expectations (vector) and SSE]
% Inputs:
%   parameters : Vector, with three values (third optional) all between 0 and 1.
%                First is learning rate second is weight of y, the third is
%                starting point if not set will be 4.5 .
%   outOther   : Vector, with outcome values for other.
%   answerSelf : Vector, answers for oneself.
%   answerOther: Vector, answers from other.

%% Function
function [sweepVec, sumSquareError] = sweepmodelself(data, param, sweepCorr, startVal)
    learnVec = zeros(length(data),1);
    sweepVec = zeros(length(data),1);
    resetVal = data(:,3);
    for i = 1:length(data)

        % For each new profile you need to reset
        if(resetVal(i) == 1)
            if(length(param) ~= 3)
                spaceVec = zeros( length(sweepCorr), 1) +  startVal;
            else
                spaceVec = zeros( length(sweepCorr), 1) +  param(3);
            end
        end

        learnVec(i) = spaceVec(data(i,10));           
        delta = data(i,5) - learnVec(i);

        spaceVec = spaceVec + param(1) * delta * sweepCorr(data(i,10));

        sweepVec(i) = learnVec(i) * param(2) + (1-param(2))*data(i,6);
    end
    sumSquareError = sum( ( data(:,4) - sweepVec ).^2 );
end