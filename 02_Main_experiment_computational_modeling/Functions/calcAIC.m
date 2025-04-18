%% AIC function
% AIC= nxLn(RSS/n)+2xK
% Where Ln is the natural log, 
%RSS is the residual sum of squares, 
%n is the number of observations that make up the sample, 
%and K is the number of parameters in your model. 

function akaikeValue = calcAIC(numTrials, numParam, SSE)
    akaikeValue =  numTrials * (log((SSE)/numTrials)) + 2 * numParam;
   
end
