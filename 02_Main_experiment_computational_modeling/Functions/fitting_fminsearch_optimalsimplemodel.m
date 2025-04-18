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
function [estimates, model, fval, exitflag, resultVec] = fitting_fminsearch_optimalsimplemodel(data, param, startVal)
%     if nargin == 3
%         parameters = param(1);
%     else
%         parameters = param([1,3]);
%     end
    parameters = param;
    %global startVal
    model = @optimalsimplemodel;
    options = optimset('Display', 'final', 'MaxIter', 10000000);
    [estimates, fval, exitflag] = fminsearch(model, parameters, options);     

    function sumSquareError = optimalsimplemodel(model_parameters)
        resultVec = zeros(length(data), 1);
        
         % Actual model
         for iModel = 1:length(data)
             
             if(data(iModel,3) == 1)
               if(length(model_parameters) ~= 2)
                   resultVec(iModel) = startVal;
               else
                   resultVec(iModel) = model_parameters(2);
               end
             end
             
             % Calculation of Prediction Error.
             delta = data(iModel,5) - resultVec(iModel);
             resultVec(iModel + 1) = resultVec(iModel) + model_parameters(1) * delta;
             
         end
         resultVec(end) = [];

        if model_parameters(1) < 0 || model_parameters(1) > 1
            sumSquareError = 10^25;
        else
            sumSquareError = sum( ( data(:,4) - resultVec ).^2 );
        end
    end

end % End function
