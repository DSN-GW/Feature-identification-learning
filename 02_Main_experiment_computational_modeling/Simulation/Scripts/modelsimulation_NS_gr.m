%% Feature identification learning 
%% Computational modeling analyses for Simulated data
%% Get actual profiles & construct parameters


addpath(genpath('./02_Main_experiment_computational_modeling/Simulation'));
addpath(genpath('./02_Main_experiment_computational_modeling/Functions'));

load ./02_Main_experiment_computational_modeling/Data/data_model_fit_all13_v1cat.mat
%dataSort=dataSort_lpr;

numPar = length(nsDataSort.data_all_cat);

%addpath(genpath('/Users/grosenblau/work_stuff_to_backup/preferences_behave/chris_matlab/revision_BioPsych/learning_Koen_HH_struct'));
%addpath(genpath('/Users/grosenblau/work_stuff_to_backup/spm12'));
%% 

%% Modelling

alphaParam = 0.2:0.05:0.8;
gammaParam = 0.2:0.05:0.8;
startValOptimal = 5.5;

repeats = length(alphaParam) * length(gammaParam);
rootRep = sqrt(repeats);
paramSimul = zeros(repeats,2);
for iAlpha = 1:rootRep
    for iGamma = 1:rootRep
        paramSimul( (iAlpha-1) * rootRep + iGamma, 1) = alphaParam(iAlpha);
        paramSimul( (iAlpha-1) * rootRep + iGamma, 2) = gammaParam(iGamma);
    end
end

%% Run the profiles for all the parameter settings
dim = 5;
% numPar = length(dataSort.delBothSortBoth);
perPart = cell(1, numPar);

for iPar = 1:numPar


  % order of the data file
       % 1-profile
       % 2-trialnr
       % 3-reset
       % 4-rating
       % 5-outcome//feedback
       % 6-self rating-c
       % 7-self rating-l
       % 8-mean rating-c
       % 9-mean rating-l
       %10-item nr (only for sweep)

    
       data = table2array(nsDataSort.data_all_cat{1,iPar}(:,[6,4,7,2,5,8:11])); % take cat or subcat sort
       data(data(:,1)==2,6) = data(data(:,1)==2,7); % own rating
       data(data(:,1)==2,8) = data(data(:,1)==2,9); % mean rating
              
       dataSweep = table2array(nsDataSort.data_all{1,iPar}(:,[6,4,7,2,5,8:11,14])); % take sweep sort
      
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
        % Pseudo linear regression
        simMat(:, 1, iSimul) = pseudolinearregressionmodel(data, paramSimul(iSimul,1:2));

        % rl model
        simMat(:, 2, iSimul) = optimalsimplemodel(data, paramSimul(iSimul, 1), startValOptimal);

        % rl self model
        simMat(:, 3, iSimul) = rlselfmodel(data, paramSimul(iSimul,1:2), startValOptimal);
        
        % Sweep model
        simMat(:, 4, iSimul) = sweepmodel(dataSweep, paramSimul(iSimul,1), sweepR, startValOptimal);

         % Sweep model self
        simMat(:, 5, iSimul) = sweepmodelself(dataSweep, paramSimul(iSimul,1:2), sweepR, startValOptimal);

%         % rl model sweep sort
%         simMat(:, 6, iSimul) = optimalsimplemodel(dataSweep, paramSimul(iSimul,1), startValOptimal);
% 
%         % rl self model sweep sort
%         simMat(:, 7, iSimul) = rlselfmodel(dataSweep, paramSimul(iSimul,1:2), startValOptimal);
%         
%         % weighted model
%         simMat(:, 8, iSimul) = rlselfmodel_weighted(dataSweep, paramSimul(iSimul,1:2), startValOptimal);

       

    end
    perPart{iPar} = simMat;
end

%% Construct noise vector and add it the the outcomes (noise vector should be size of iterations wanted)
%% Construct noise vector and add it to the outcomes (noise vector should be size of iterations wanted)
numIt = 200;
noisePart = cell(numIt,numPar);
for it = 1:numIt
    % Add noise to outcome matrix
    for iAdd = 1:numPar
        noiseMat = randn(size(perPart{iAdd}));
        noisePart{it,iAdd} = perPart{iAdd} + noiseMat;
        noisePart{it,iAdd}(noisePart{it,iAdd} < 0) = 0; %% this should be bound better...
    end
end

%% Run fminsearch over all these outcomes with noise
% Need to exchange the columns with those I just created.
% param    = [0.5, 0.5, 5.5];
% startVal = 5.5;
param       = [0.5, 0.5, 5.5];
paramSimple = [0.5, 5.5];
paramRecover = repmat([paramSimul(:,1), zeros(repeats,1), paramSimul(:,2), zeros(repeats,1)],1,1,dim,numPar);

%cat or subcat: change this
which_cat = 'cat' %'sub_cat' %'cat'
data_nsSort=nsDataSort.(subsref(fieldnames(nsDataSort),substruct('{}',{3})));
bicMat   = zeros(dim);

for iIt = 1:numIt
    for iParam = 1:repeats
        for iModel = 1:dim
            bicDim = zeros(1, dim);
            for iPar = 1:numPar

            data = table2array(data_nsSort{1,iPar}(:,[6,4,7,2,5,8:11])); % take cat or subcat sort
            data(data(:,1)==2,6) = data(data(:,1)==2,7); % own rating
            data(data(:,1)==2,8) = data(data(:,1)==2,9); % mean rating
            data(:,4) = noisePart{iIt,iPar}(:,iModel,iParam);
            
            dataSweep = table2array(nsDataSort.data_all{1,iPar}(:,[6,4,7,2,5,8:11,14])); % take sweep sort
            dataSweep(dataSweep(:,1)==2,6) = dataSweep(dataSweep(:,1)==2,7); % own rating
            dataSweep(dataSweep(:,1)==2,8) = dataSweep(dataSweep(:,1)==2,9); % mean rating
            dataSweep(:,4) = noisePart{iIt,iPar}(:,iModel,iParam);

            % convert sweep matrix number 
            for i_tem=1:length(dataSweep)
                if sum(dataSweep(i_tem,10)==convert(:,1))==1
                dataSweep(i_tem,10)=convert(find(dataSweep(i_tem,10)==convert(:,1)),2);
                end
           end

            trials = length(data);    

           
            % ### pseudo linear regression model ###
            [estimates, ~, fval, ~, ~] = fitting_fminsearch_pseudolinearregression(data, paramSimple);
            paramRecover(iParam,2,1,iPar) = estimates(1);
            if length(estimates) == 3
                paramRecover(iParam,4,1,iPar) = estimates(2);
            end
             bicDim(1) = bicDim(1) + calcBIC(trials, 2, fval);
            
            
            % ### rl model ###
            [estimates, ~, fval, ~, ~] = fitting_fminsearch_optimalsimplemodel(data, paramSimple);
            paramRecover(iParam,2,2,iPar) = estimates(1);
            if length(estimates) == 3
                paramRecover(iParam,4,2,iPar) = estimates(2);
            end
            bicDim(2) = bicDim(2) + calcBIC(trials, 2, fval);

           
            % ### rl self model ###
            [estimates, ~, fval, ~, ~] = fitting_fminsearch_rlselfmodel(data, param);
            paramRecover(iParam,2,3,iPar) = estimates(1);
            if length(estimates) == 3
                paramRecover(iParam,4,3,iPar) = estimates(2);
            end
            bicDim(3) = bicDim(3) + calcBIC(trials, 3, fval);
            
            
            % ### sweep model ###
            [estimates, ~, fval, ~, ~,~] = fitting_fminsearch_sweepmodel(dataSweep, paramSimple, sweepR);
            paramRecover(iParam,2,4,iPar) = estimates(1);
            if length(estimates) == 3
                paramRecover(iParam,4,4,iPar) = estimates(2);
            end
            bicDim(4) = bicDim(4) + calcBIC(trials, 2, fval);


            % ### sweep self model ###
            [estimates, ~, fval, ~, ~,~] = fitting_fminsearch_sweepmodelself(dataSweep, param, sweepR);
            paramRecover(iParam,2,5,iPar) = estimates(1);
            if length(estimates) == 3
                paramRecover(iParam,4,5,iPar) = estimates(2);
            end
            bicDim(5) = bicDim(5) + calcBIC(trials, 3, fval);


%             dataSweep(:,4) = noisePart{iIt,iPar}(:,6,iParam);
%             % ### rl model sweep sort ###
%             [estimates, ~, fval, ~, ~] = fitting_fminsearch_optimalsimplemodel(dataSweep, param);
%             paramRecover(iParam,2,6,iPar) = estimates(1);
%             if length(estimates) == 3
%                 paramRecover(iParam,4,6,iPar) = estimates(2);
%             end
%             bicDim(6) = bicDim(6) + calcBIC(trials, 2, fval);
% 
%              dataSweep(:,4) = noisePart{iIt,iPar}(:,7,iParam);
%             % ### rl self model sweep sort ###
%             [estimates, ~, fval, ~, ~] = fitting_fminsearch_rlselfmodel(dataSweep, param);
%             paramRecover(iParam,2,7,iPar) = estimates(1);
%             if length(estimates) == 3
%                 paramRecover(iParam,4,7,iPar) = estimates(2);
%             end
%             bicDim(7) = bicDim(7) + calcBIC(trials, 3, fval);
% 
%             dataSweep(:,4) = noisePart{iIt,iPar}(:,8,iParam);
%             % ### rl self model weighted ###
%             [estimates, ~, fval, ~, ~,~] = fitting_fminsearch_rlselfmodel_weighted(dataSweep, param);
%             paramRecover(iParam,2,8,iPar) = estimates(1);
%             if length(estimates) == 3
%                 paramRecover(iParam,4,8,iPar) = estimates(2);
%             end
%             bicDim(8) = bicDim(8) + calcBIC(trials, 3, fval);

            end
            bicMat(iModel, :) = bicMat(iModel, :) + (bicDim == min(bicDim));
        end
    end
end

                
%          


%%
save ('simulation/modelsim_simpar.mat');

%% Plot parameter values (scatterplot) & plot Confusion Matrix (imagesc) (Maybe in grey values? grey to black)
figure
imagesc(bicMat)
colormap('gray')

%% check corrs of model predictions
for iPar = 1:numPar
c{iPar} = corrcoef(noisePart{1,1}(:,[1 2:5])); 
end

%% Plot parameter values (scatterplot)
recover = [];
for iAdd = 1:numPar
    recover = [recover; paramRecover(:,:,:,iPar)];
end

recover_mean = mean(paramRecover,4);


 for iPlot = [4 5]
    figure
    scatter(paramRecover(:,1,iPlot),paramRecover(:,2,iPlot))
    title(['Model: ',num2str(iPlot), ' First parameter']);
    xlabel('Input')
    ylabel('Recovered')
    figure
    scatter(paramRecover(:,3,iPlot),paramRecover(:,4,iPlot))
     title(['Model: ',num2str(iPlot), ' Second parameter']);
end

%% check corrs

corr(recover_mean(:,1,4),recover_mean(:,2,4))
corr(recover_mean(:,1,5),recover_mean(:,2,5))
corr(recover_mean(:,3,5),recover_mean(:,4,5))

    figure
    scatter(recover_mean(:,1,4),recover_mean(:,2,4))
    title(['Model: ',num2str(5), ' First parameter']);
    xlabel('Input')
    ylabel('Recovered')
%     figure
%     scatter(paramRecover(:,3,iPlot),paramRecover(:,4,iPlot))
%      title(['Model: ',num2str(iPlot), ' Second parameter']);



figure
for iP = 1:dim
    subplot(2,4,iP)
    scatter(recover_mean(:,1,iP),recover_mean(:,2,iP)) 

%     if iP==3
%     subplot(2,4,iP)
%     scatter(recover_mean(:,1,iP),recover_mean(:,2,iP)) 
%     subplot(2,4,iP+1)
%     scatter(recover_mean(:,3,iP),recover_mean(:,4,iP))
%   
%     
%     elseif iP==4
%     subplot(2,4,iP+1)
%     scatter(recover_mean(:,1,4),recover_mean(:,2,4))  
% 
% 
%     elseif iP==5
%     subplot(2,4,iP+1)
%     scatter(recover_mean(:,1,5),recover_mean(:,2,5))    
%     subplot(2,4,iP+2)
%     scatter(recover_mean(:,3,5),recover_mean(:,4,5))
%     end

    xlabel('Input')
    ylabel('Recovered')
    title(['Model: ', num2str(iP)])
end

figure
for iP = 1:dim
    subplot(2,3,iP)
    scatter(recover_mean(:,3,iP),recover_mean(:,4,iP))
    xlabel('Input')
    ylabel('Recovered')
    title(['Model: ', num2str(iP)])
end


for iP = 1:dim
    c{iP} = corrcoef(recover_mean(:,1,iP),recover_mean(:,2,iP))
end

for iP = [3 5]
    c2{iP} = corrcoef(recover_mean(:,1,iP),recover_mean(:,2,iP))
end

MeanC = mean(cat(3, c{:}), 3);  
figure
imagesc(tril(MeanC))
cb = colorbar; 
colormap default;
ylabel(cb,'correlation coefficient')


%% BIC function
function bayesianValue = calcBIC(numTrials, numParam, SSE)
    bayesianValue = ( numTrials * log(SSE/numTrials) ) + ( numParam * log(numTrials) );
end
