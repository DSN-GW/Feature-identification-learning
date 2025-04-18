%% rlselfmodel
% combination of RL & self-rating according to formula by Y. Niv: 
% rating = a*Vt + (1-a)own

%% Function
function [learnAndSelfVec, sumSquareError] = rlselfmodel(data, parameters, startVal)
    learnVec        = zeros(length(data),1);
    learnAndSelfVec = zeros(length(data),1);
    %global startVal

    for i = 1:length(data)
         % This is the reset for a new profile, resets are column 5
         if(data(i,3) == 1)

             if(length(parameters) ~= 3); learnVec(i) = startVal
             elseif(length(parameters) == 3); learnVec(i) = parameters(3); end

         end

        delta = data(i,5) - learnVec(i);
        learnVec(i+1) = learnVec(i) + parameters(1) * delta;
        learnAndSelfVec(i) = parameters(2) * learnVec(i) + (1 - parameters(2)) * data(i,6);
    end
    
    % Since we use fminsearch (No bound)
    if parameters(1) < 0 || parameters(1) > 1
        sumSquareError = 10^25;
    elseif parameters(2) < 0 || parameters(2) > 1
        sumSquareError = 10^25;
    else
        sumSquareError = sum( ( data(:,4) - learnAndSelfVec ).^2 );
    end
end % End rlselfmodel