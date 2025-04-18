
%% Feature identification learning 
%% plot results from computational modeling
%%
% #########################################################################

%% cat or subcat: change this to change to category resets
which_cat = 'sub_cat' %'sub_cat' %'cat'

load(fullfile('./02_Main_experiment_computational_modeling/Data/data',strcat('data_model_fit_all13_v1',which_cat, '.mat')));
addpath('./02_Main_experiment_computational_modeling/Functions/spm12')
cd ('./02_Main_experiment_computational_modeling/Script/plotting')
%% plot data for main model types
%% make sure that SPM is in your path

%% plots the main model space
figTitle = 'data 5 models';
[BMS] = PXP_random_plot_v1( -bicMatrix(:,[1,3,4,6,7]), 'PXP random plot' );
newYticklabel = {'Fine Granularity & Reference Point RP','Fine Granularity', 'Coarse Granularity & Reference Point RP','Coarse Granularity' ,'No learning' };
set(gca,'YtickLabel',newYticklabel);
BIC_fixed_plot_v1( bicMin(:,[1,3,4,6,7]), ['BIC scores: ',figTitle])
newYticklabel = {'Fine Granularity & Reference Point RP','Fine Granularity', 'Coarse Granularity & Reference Point RP','Coarse Granularity' ,'No learning' };
set(gca,'YtickLabel',newYticklabel);

%% parameter plots

Like2=bicMatrix(:,[1,3,4,6,7])'
numPar = length(nsDataSort.data_all);
% plot frequencies of participants using a certain model
for i_par=1:numPar
[row(i_par), column(i_par)] = find(Like2(:,i_par) == min(Like2(:,i_par)));
end

figure
h1=axes
[counts, centers] = histcounts(row);
barh(counts)
set(h1, 'Ydir', 'reverse')

%% plot parameters
for iPar=1:numPar
estimates(iPar,1) = paramEstimates( 6 ).model{1,iPar}(1);
estimates(iPar,2) = paramEstimates( 6 ).model{1,iPar}(2);
end

figure
subplot(2,1,1)
hist( estimates(:,1) )
xlim([0 1])

subplot(2,1,2)
hist( estimates(:,2)  )
xlim([1 10])

figure 
boxplot([estimates(:,1),estimates(:,2)])
newXticklabel = {'alpha', 'gamma' };
set(gca,'XtickLabel',newXticklabel);

%% plot similarity matrix
B = rescale(sweepR)
figure
imagesc(B)
clims = [0.9 1];
imagesc(sweepR)

