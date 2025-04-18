%% BIC function
function bayesianValue = calcBIC(numTrials, numParam, SSE)
    bayesianValue = ( numTrials * log(SSE/numTrials) ) + ( numParam * log(numTrials) );
end
