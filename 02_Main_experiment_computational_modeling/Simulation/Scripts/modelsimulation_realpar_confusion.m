%% Feature identification learning 
%% Computational modeling analyses for Simulated data
%% Get actual profiles & construct parameters

load ./Nonsocial_learning_task/modeling/data/data_model_fit_all13_v1cat.mat

addpath(genpath('./Nonsocial_learning_task/modeling/simulation'));
addpath(genpath('./Nonsocial_learning_task/modeling/fitting_functions'));
%% 

%% Modelling


alphaParam = 1;
gammaParam = 1;
startValOptimal = 5.5;
repeats = 1;

%% Run the profiles for all the parameter settings
dim = 5;
 numPar = length(nsDataSort.data_all_cat);
perPart = cell(1, numPar);
for iPar = 1:numPar
    
       data = table2array(nsDataSort.data_all_cat{1,iPar}(:,[6,4,7,2,5,8:11])); % take cat or subcat sort
       data(data(:,1)==2,6) = data(data(:,1)==2,7); % own rating
       data(data(:,1)==2,8) = data(data(:,1)==2,9); % mean rating
              
       dataSweep = table2array(nsDataSort.data_all{1,iPar}(:,[6,4,7,2,5,8:11,14])); % take sweep sor
       dataSweep(dataSweep(:,1)==2,6) = dataSweep(dataSweep(:,1)==2,7); % own rating
       dataSweep(dataSweep(:,1)==2,8) = dataSweep(dataSweep(:,1)==2,9); % mean rating

     % convert sweep matrix number 
       for i_tem=1:length(dataSweep)
           if sum(dataSweep(i_tem,10)==convert(:,1))==1
               
               dataSweep(i_tem,10)=convert(find(dataSweep(i_tem,10)==convert(:,1)),2);
               
           end
       end
    

    trials = length(data);
    simMat = zeros(trials, dim, length(alphaParam)*length(gammaParam));
    for iSimul = 1:repeats
        
        %regression
        simMat(:, 1, iSimul) = pseudolinearregressionmodel(data, paramEstimates(1).model{1,iPar});

        % rl model
        simMat(:, 2, iSimul) = optimalsimplemodel(data, paramEstimates(3).model{1,iPar});

        % rl self model
        simMat(:, 3, iSimul) = rlselfmodel(data, paramEstimates(4).model{1,iPar});
        
        % Sweep model
        simMat(:, 4, iSimul) = sweepmodel(dataSweep, paramEstimates(6).model{1,iPar}, sweepR);

         % Sweep model self
        simMat(:, 5, iSimul) = sweepmodelself(dataSweep, paramEstimates(7).model{1,iPar}, sweepR);

        % rl model sweep sort
%         simMat(:, 6, iSimul) = optimalsimplemodel(dataSweep, paramEstimates(9).model{1,iPar});
% 
%         % rl self model sweep sort
%         simMat(:, 7, iSimul) = rlselfmodel(dataSweep, paramEstimates(10).model{1,iPar});
% 
%         % weighted model
%         simMat(:, 8, iSimul) = rlselfmodel_weighted(dataSweep, paramEstimates(12).model{1,iPar});

        
    end
    perPart{iPar} = simMat;
end

%% Construct noise vector and add it the the outcomes (noise vector should be size of iterations wanted)
numIt = 200;
noisePart = cell(numIt,numPar);
for it = 1:numIt
    % Add noise to outcome matrix
    for iAdd = 1:numPar
        noiseMat = randn(size(perPart{iAdd}));
        noisePart{it,iAdd} = perPart{iAdd} + noiseMat;
        noisePart{it,iAdd}(noisePart{it,iAdd} < 1) = 1;
        noisePart{it,iAdd}(noisePart{it,iAdd} > 10) = 10;
    end
end

%% Run fminsearch over all these outcomes with noise
% Need to exchange the columns with those I just created.
% param    = [0.5, 0.5, 5.5];
% startVal = 5.5;

param       = [0.5, 0.5, 5.5];
paramSimple = [0.5, 5.5];
bicMat   = zeros(dim);
%bicMat   = zeros(7);
totParam = [];

for iIt = 1:numIt
    for iParam = 1:repeats
        for iModel = 1:dim
            bicDim = zeros(1, dim);
            for iPar = 1:numPar
                
           data = table2array(nsDataSort.data_all_cat{1,iPar}(:,[6,4,7,2,5,8:11])); % take cat or subcat sort
           data(data(:,1)==2,6) = data(data(:,1)==2,7); % own rating
           data(data(:,1)==2,8) = data(data(:,1)==2,9); % mean rating
           data(:,4) = noisePart{iIt,iPar}(:,iModel,iParam);


           dataSweep = table2array(nsDataSort.data_all{1,iPar}(:,[6,4,7,2,5,8:11,14])); % take sweep sor
           dataSweep(dataSweep(:,1)==2,6) = dataSweep(dataSweep(:,1)==2,7); % own rating
           dataSweep(dataSweep(:,1)==2,8) = dataSweep(dataSweep(:,1)==2,9); % mean rating
           dataSweep(:,4) = noisePart{iIt,iPar}(:,iModel,iParam);

          % convert sweep matrix number 
            for i_tem=1:length(dataSweep)
            if sum(dataSweep(i_tem,10)==convert(:,1))==1
               
               dataSweep(i_tem,10)=convert(find(dataSweep(i_tem,10)==convert(:,1)),2);
               
            end
            end


%                 data = dataSort.delBothSortBoth{iPar};
                trials = length(data);

                % ### pseudo linear regression model ###     
                [param_Esti(1,iModel).model{iIt,iPar}, ~, fval, ~, ~] = fitting_fminsearch_pseudolinearregression(data, paramSimple);
                bicDim(1) = bicDim(1) + calcBIC(trials, 2, fval);

               % ### rl model ###
                [param_Esti(2,iModel).model{iIt,iPar}, ~, fval, ~, ~] = fitting_fminsearch_optimalsimplemodel(data, paramSimple);
                bicDim(2) = bicDim(2) + calcBIC(trials, 2, fval);

                % ### rl self model ###
                [param_Esti(3,iModel).model{iIt,iPar}, ~, fval, ~, ~] = fitting_fminsearch_rlselfmodel(data, param);          
                bicDim(3) = bicDim(3) + calcBIC(trials, 3, fval);
                
                % ### rl sweep ###
                [param_Esti(4,iModel).model{iIt,iPar}, ~, fval, ~, ~] = fitting_fminsearch_sweepmodel(dataSweep, paramSimple, sweepR);
                bicDim(4) = bicDim(4) + calcBIC(trials, 2, fval);
                
                % ### sweep comb ###
            
                [param_Esti(5,iModel).model{iIt,iPar}, ~, fval, ~, ~] = fitting_fminsearch_sweepmodelself(dataSweep, param, sweepR);
                bicDim(5) = bicDim(5) + calcBIC(trials, 3, fval);

%                 % ### RL sweep sort model ###
%                 [param_Esti(6,iModel).model{iIt,iPar}, ~, fval, ~, ~] = fitting_fminsearch_optimalsimplemodel(dataSweep, paramSimple);
%                 bicDim(6) = bicDim(6) + calcBIC(trials, 2, fval);
% 
%                 % ### RL comb sweep sort ###
%                 [param_Esti(7,iModel).model{iIt,iPar}, ~, fval, ~, ~, ~] = fitting_fminsearch_rlselfmodel(dataSweep, param);
%                 bicDim(7) = bicDim(7) + calcBIC(trials, 3, fval);
% 
%                 % ### RL weight ###
%                 [param_Esti(8,iModel).model{iIt,iPar}, ~, fval, ~, ~, ~] = fitting_fminsearch_rlselfmodel_weighted(dataSweep, param);
%                 bicDim(8) = bicDim(8) + calcBIC(trials, 3, fval);
%             
            
            end
            bicMat(iModel, :) = bicMat(iModel, :) + (bicDim == min(bicDim));
        end
    end
end

save ('./Nonsocial_learning_task/modeling/simulation/modelsim_realpar.mat');

%% Plot parameter values (scatterplot) & plot Confusion Matrix (imagesc) (Maybe in grey values? grey to black)
figure
imagesc(bicMat)
colormap('gray')



recover = [];
for iAdd = 1:numPar
    recover = [recover; param_Esti(:,:,:,iPar)];
end
paramRecover(iParam,2,dim,iPar)

 for iPlot = 1:dim
    figure
    scatter(paramRecover(:,1,iPlot),paramRecover(:,2,iPlot))
    title(['Model: ',num2str(iPlot), ' First parameter']);
    xlabel('Input')
    ylabel('Recovered')
    figure
    scatter(paramRecover(:,3,iPlot),paramRecover(:,4,iPlot))
     title(['Model: ',num2str(iPlot), ' Second parameter']);
end




%% check corrs of model predictions
for iPar = 1:numPar
c{iPar} = corrcoef(noisePart{1,1}(:,[1 2:7])); 
end


MeanC = mean(cat(3, c{:}), 3);  % Or: mean(C{iC}(:)) ?
figure
imagesc(tril(MeanC))
cb = colorbar; 
colormap default;
ylabel(cb,'correlation coefficient')




%% check par corr

param_recover{8} = {};

for imod_sim =1:7 %sim data
for imod=1:7 % running model
    for iPar = 1:numPar
        for iIt = 1:numIt
           inst(iIt,:) = param_Esti(imod,imod_sim).model{iIt,iPar}(1:2);
        end    
        param_recover{imod,imod_sim}(iPar,:) = mean(int);
        param_recover{8,imod_sim}(iPar,:) = paramEstimates(1).model{1,iPar}(1:2);
        param_recover{9,imod_sim}(iPar,:) = paramEstimates(2).model{1,iPar}(1:2);
        param_recover{10,imod_sim}(iPar,:) = paramEstimates(3).model{1,iPar}(1:2);
        param_recover{11,imod_sim}(iPar,:) = paramEstimates(4).model{1,iPar}(1:2);
        param_recover{12,imod_sim}(iPar,:) = paramEstimates(5).model{1,iPar}(1:2);
        param_recover{13,imod_sim}(iPar,:) = paramEstimates(6).model{1,iPar}(1:2);
        param_recover{14,imod_sim}(iPar,:) = paramEstimates(7).model{1,iPar}(1:2);
    end
    clear int
    end
end


mod1_recover = ( corrcoef([param_recover{1,1}(:,1)...
          param_recover{1,2}(:,1)... 
          param_recover{1,3}(:,1)...
          param_recover{1,4}(:,1)...
          param_recover{1,5}(:,1)... 
          param_recover{1,6}(:,1)...
          param_recover{1,7}(:,1)...
          param_recover{8,1}(:,1)]));

mod2_recover = ( corrcoef([param_recover{2,1}(:,1)...
          param_recover{2,2}(:,1)... 
          param_recover{2,3}(:,1)...
          param_recover{2,4}(:,1)...
          param_recover{2,5}(:,1)... 
          param_recover{2,6}(:,1)...
          param_recover{2,7}(:,1)...
          param_recover{9,1}(:,1)]));
      
mod3_recover = ( corrcoef([param_recover{3,1}(:,1)...
          param_recover{3,2}(:,1)... 
          param_recover{3,3}(:,1)...
          param_recover{3,4}(:,1)...
          param_recover{3,5}(:,1)... 
          param_recover{3,6}(:,1)...
          param_recover{3,7}(:,1)...
          param_recover{10,1}(:,1)]));
      
 mod4_recover = ( corrcoef([param_recover{4,1}(:,1)...
          param_recover{4,2}(:,1)... 
          param_recover{4,3}(:,1)...
          param_recover{4,4}(:,1)...
          param_recover{4,5}(:,1)... 
          param_recover{4,6}(:,1)...
          param_recover{4,7}(:,1)...
          param_recover{11,1}(:,1)]));
      
 mod5_recover = ( corrcoef([param_recover{5,1}(:,1)...
          param_recover{5,2}(:,1)... 
          param_recover{5,3}(:,1)...
          param_recover{5,4}(:,1)...
          param_recover{5,5}(:,1)... 
          param_recover{5,6}(:,1)...
          param_recover{5,7}(:,1)...
          param_recover{12,1}(:,1)]));

 mod6_recover = ( corrcoef([param_recover{6,1}(:,1)...
          param_recover{6,2}(:,1)... 
          param_recover{6,3}(:,1)...
          param_recover{6,4}(:,1)...
          param_recover{6,5}(:,1)... 
          param_recover{6,6}(:,1)...
          param_recover{6,7}(:,1)...
          param_recover{13,1}(:,1)]));

 mod7_recover = ( corrcoef([param_recover{7,1}(:,1)...
          param_recover{7,2}(:,1)... 
          param_recover{7,3}(:,1)...
          param_recover{7,4}(:,1)...
          param_recover{7,5}(:,1)... 
          param_recover{7,6}(:,1)...
          param_recover{7,7}(:,1)...
          param_recover{14,1}(:,1)]));
figure
imagesc(abs([mod1_recover(8,1:7); mod2_recover(8,1:7);...
        mod3_recover(8,1:7); mod4_recover(8,1:7);...
        mod5_recover(8,1:7); mod6_recover(8,1:7);...
        mod7_recover(8,1:7)]))
colormap default

save('./Nonsocial_learning_task/modeling/simulation/modelsim_.mat' )


%% BIC function
function bayesianValue = calcBIC(numTrials, numParam, SSE)
    bayesianValue = ( numTrials * log(SSE/numTrials) ) + ( numParam * log(numTrials) );
end
