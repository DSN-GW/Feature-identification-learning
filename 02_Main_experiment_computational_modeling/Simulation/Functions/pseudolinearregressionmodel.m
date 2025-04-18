%% fitting_fminsearch_pseudolinearregression
% Should work like a linear regression

function [resultVec,sumSquareError] = pseudolinearregressionmodel(data, model_parameters)
    resultVec = model_parameters(1) * data(:,6) + model_parameters(2);
    sumSquareError = sum( (data(:,4) - resultVec).^2 );
end