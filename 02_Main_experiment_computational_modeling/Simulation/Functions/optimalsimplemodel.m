%% fitting_fminsearch_optimalsimplemodel
% Simplest RL model independent from participant response
%
% Output is a vector with predictions
% Inputs:
%   learnRate: Single value, between 0 and 1.
%   outcome  : Vector, with outcome values.
%   actual   : Vector, with "actual" outcome values.
%   startVal : (optional)Single value, if not set will be 0.

%% Function
function [resultVec, sumSquareError] = optimalsimplemodel(data, parameters, startVal)
    resultVec = zeros(length(data), 1);

     % Actual model
     for iModel = 1:length(data)

         % This is the reset for a new profile, resets are column 21
         if(data(iModel,3) == 1)
             if(length(parameters) ~= 2); resultVec(iModel) = startVal
             %resultVec(iModel) = startVal;
             elseif(length(parameters) == 2); resultVec(iModel) = parameters(2); end
         end

         % Calculation of Prediction Error.
         delta = data(iModel,5) - resultVec(iModel);
         resultVec(iModel + 1) = resultVec(iModel) + parameters(1) * delta;

     end
     resultVec(end) = [];

     sumSquareError = sum( ( data(:,4) - resultVec ).^2 );
end % End function