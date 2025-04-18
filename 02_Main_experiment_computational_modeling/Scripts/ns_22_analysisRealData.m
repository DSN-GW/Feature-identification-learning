%% Feature identification learning 
%% Computational modeling analyses for Real Data sets 
% #########################################################################

clear all;
close all;

%% adding path 
addpath(genpath('./Nonsocial_learning_task/modeling'));
%% add SPM to your path
addpath(genpath('./spm12'));

%% Load data 
cd ('./');
load ('./02_Main_experiment_computational_modeling/Data/nsDataSortModel.mat','nsDataSort' );
load ('./02_Main_experiment_computational_modeling/Data/item_distances_exp1.mat');
convert=readtable ('./02_Main_experiment_computational_modeling/Data/conversions.csv');
convert=table2array(convert);
numPar = length(nsDataSort.data_all);
sweepR = exp1.all_transform;
sweepR = rescale(sweepR);

%% Modelling structures
numModels      = 8;
startVal       = 5.5; % Used as a fixed value. If none set input as free parameter
param          = [0.5, 0.5, 5.5];
%results        = struct('estimates', repmat(struct('model', {cell(1,numPar)}),numModels, 1 ) );
paramEstimates = struct(repmat(struct('model', {cell(1,numPar)}),numModels, 1 ) );
regressCheck   = zeros(numPar, 2);
bicMatrix   = zeros(numPar, numModels);
aicMatrix   = zeros(numPar, numModels);
bicSummed = zeros(1,numModels);

%% cat or subcat: change this to change to category resets
which_cat = 'sub_cat' %'sub_cat' %'cat'
data_nsSort=nsDataSort.(subsref(fieldnames(nsDataSort),substruct('{}',{1})))


for iModels = 1:numModels

    
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

    
       data = table2array(data_nsSort{1,iPar}(:,[6,4,7,2,5,8:11])); % take cat or subcat sort
       data(data(:,1)==2,6) = data(data(:,1)==2,7); % own rating
       data(data(:,1)==2,8) = data(data(:,1)==2,9); % mean rating
      
       %% data using mean and not self
       dataMean = data;  
       dataMean(:,6) = data(:,8);  
%         
       dataSweep = table2array(nsDataSort.data_all{1,iPar}(:,[6,4,7,2,5,8:11,14])); % take sweep sort
       dataSweep(dataSweep(:,1)==2,6) = dataSweep(dataSweep(:,1)==2,7); % own rating
       dataSweep(dataSweep(:,1)==2,8) = dataSweep(dataSweep(:,1)==2,9); % mean rating
     
       
       % convert sweep matrix number 
       for i_tem=1:length(dataSweep)
           if sum(dataSweep(i_tem,10)==convert(:,1))==1
               
               dataSweep(i_tem,10)=convert(find(dataSweep(i_tem,10)==convert(:,1)),2);
               
           end
       end

       dataSweepMean = dataSweep;  
       dataSweepMean(:,6) = dataSweep(:,8);
           
%         dataSweepMean = dataSweep;
%         dataSweep(:,16) = dataSweep(:,4);
%         newYticklabel = {'Fine Granularity & Reference Point RP','Fine Granularity', 'Coarse Granularity & Reference Point RP','Coarse Granularity' ,'No learning' };
                    
        if iModels == 1 % Uses both sorts and both deletes.
            [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}]   = fitting_fminsearch_pseudolinearregression(data, param([1 3]));
            bicMatrix(iPar,1) = calcBIC(length(data), 2, fval);
            aicMatrix(iPar,1) = calcAIC(length(data), 2, fval);
            
        elseif iModels == 2 % Uses both sorts and both deletes.
             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}]   = fitting_fminsearch_pseudolinearregression(dataMean, param([1 3]));
            bicMatrix(iPar,2) = calcBIC(length(dataMean), 2, fval);
            aicMatrix(iPar,2) = calcAIC(length(dataMean), 2, fval);

        elseif iModels == 3 % Uses both sorts and other deletes.
             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}]   = fitting_fminsearch_optimalsimplemodel(data, param([1 3]));
            bicMatrix(iPar,3) = calcBIC(length(data), 2, fval);
            aicMatrix(iPar,3) = calcAIC(length(data), 2, fval);

        elseif iModels == 4 % Uses both sorts and both deletes.
             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_rlselfmodel(data, param);
            bicMatrix(iPar,4) = calcBIC(length(data), 3, fval);
            aicMatrix(iPar,4) = calcAIC(length(data), 3, fval);
            % ### Same model but with mean ratings
        elseif iModels == 5 % Uses both sorts and both deletes.
            [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_rlselfmodel(dataMean, param);
            bicMatrix(iPar,5) = calcBIC(length(dataMean), 3, fval);
            aicMatrix(iPar,5) = calcAIC(length(dataMean), 3, fval);

        elseif iModels == 6 % Uses prof sort and both deletes.
            [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_sweepmodel(dataSweep, param([1 3]), sweepR);
            bicMatrix(iPar,6) = calcBIC(length(dataSweep), 2, fval);
            aicMatrix(iPar,6) = calcAIC(length(dataSweep), 2, fval);

        elseif iModels == 7 % Uses prof sort and both deletes.
            [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_sweepmodelself(dataSweep, param, sweepR);
            bicMatrix(iPar,7) = calcBIC(length(dataSweep), 3, fval);
            aicMatrix(iPar,7) = calcAIC(length(dataSweep), 3, fval);
            
        elseif iModels == 8 % Uses prof sort and both deletes.
            [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_sweepmodelself(dataSweepMean, param, sweepR);
            bicMatrix(iPar,8) = calcBIC(length(dataSweep), 3, fval);
            aicMatrix(iPar,8) = calcAIC(length(dataSweep), 3, fval);
        

            
 %% these models were initially introduced to check generalizations across the entire run (resets at run level and not cat or subcat levels)
 %% produced worse fit than original model space

    %         elseif iModels == 9 % Uses both sorts and other deletes.
%             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}]   = fitting_fminsearch_optimalsimplemodel(dataSweep, param([1 3 ]));
%             bicMatrix(iPar,9) = calcBIC(length(data), 2, fval);
%             aicMatrix(iPar,9) = calcAIC(length(data), 2, fval);
%         
%         elseif iModels == 10 % Uses both sorts and both deletes.
%             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_rlselfmodel(dataSweep, param);
%             bicMatrix(iPar,10) = calcBIC(length(data), 3, fval);
%             aicMatrix(iPar,10) = calcAIC(length(data), 3, fval);
%             % ### Same model but with mean ratings
%         
%         elseif iModels == 11 % Uses both sorts and both deletes.
%             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_rlselfmodel(dataSweepMean, param);
%             bicMatrix(iPar,11) = calcBIC(length(data), 3, fval);
%             aicMatrix(iPar,11) = calcAIC(length(data), 3, fval);
%             % ### Same model but with mean ratings
% 
%         
%         elseif iModels == 12 % Uses both sorts and both deletes.
%             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_rlselfmodel_weighted(dataSweep, param);
%             bicMatrix(iPar,12) = calcBIC(length(data), 3, fval);
%             aicMatrix(iPar,12) = calcAIC(length(data), 3, fval);
%             % ### Same model but with mean ratings
% 
%         elseif iModels == 13 % Uses both sorts and both deletes.
%             [paramEstimates(iModels).model{iPar}, model,fval,exitflag(iPar,iModels),resultV(iModels).model{iPar}] = fitting_fminsearch_rlselfmodel_weighted(dataSweepMean, param);
%             bicMatrix(iPar,13) = calcBIC(length(data), 3, fval);
%             aicMatrix(iPar,13) = calcAIC(length(data), 3, fval);
%             % ### Same model but with mean ratings


        end % End model check

    end % End iParticipants
    
end % End iModels

bicMin = sum(bicMatrix(:,:));
aicMin = sum(aicMatrix(:,:));

%% plot
% figure
% bar(bicMin)

%% aic weights
aic_delta = aicMin - min(aicMin);
[L,W]=weighted_aic(aic_delta);


%VBA package

[posterior,out] = VBA_groupBMC(Like2)

clearvars -except bestMod bicMatrix nsDataSort data_nsSort paramEstimates bicMin aicMin which_cat sweepR convert
save(fullfile('data',strcat('data_model_fit_all13_v1',which_cat, '.mat')));
