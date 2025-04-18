%% fitting_fminsearch_rlselfmodel
% combination of RL & self-rating according to formula by Y. Niv: 
% rating = a*Vt + (1-a)own
%
% Output (two most important):
%   estimates : The estimates parameter values
%   fval      : The SSE, used again for calculating the BIC
%
% Inputs:
%   data    : matrix with all neccesary data
%   param   : starting parameter values, with(out) a third starting value
%   startVal: optional starting value if you did not add a third parameter

%% Function
function [estimates, model, fval, exitflag, learnVec, learnAndSelfVec] = fitting_fminsearch_rlselfmodel(data, parameters,startVal)
    % If a startValue is set use that as a fixed value.
    % Need to set here because parameters will be varied.
%     if nargin ~= 3
%         parameters = param(1:2);
%     else
%         parameters = param;
%     end
   % global startVal
    model = @rlselfmodel;
    options = optimset('Display', 'final', 'MaxIter', 10000000);
    [estimates,fval,exitflag] = fminsearch(model, parameters, options);
   
    % Actual Model
    function [sumSquareError] = rlselfmodel(mod_param)
        learnVec        = zeros(length(data),1);    
        learnAndSelfVec = zeros(length(data),1);
        for i = 1:length(data)
              
            % First do a check for resets.
            if(data(i,3) == 1)
               if(length(mod_param) ~= 3)
                  learnVec(i) = startVal;
               else
                   learnVec(i) = mod_param(3);
               end
            end
            
            delta = data(i,5) - learnVec(i);
            learnVec(i+1) = learnVec(i) + mod_param(1) * delta;
            learnAndSelfVec(i) = mod_param(2) * learnVec(i) + (1 - mod_param(2)) * data(i,6);
        end
        
        % Since fminsearch has no bounds, do it ourselves.
        if(mod_param(1) < 0 || mod_param(1) > 1)
            sumSquareError = 10e10;
        elseif(mod_param(2) < 0 || mod_param(2) > 1)
            sumSquareError = 10e10;
        else
            sumSquareError = sum( ( data(:,4) - learnAndSelfVec ).^2 );
        end
        
    end % End rlselfmodel

end % End fminsearch