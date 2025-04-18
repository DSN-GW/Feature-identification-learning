%% Feature identification learning 
%% This script organizes data for modeling

%% column info

%  1 - subject number
%  2 - profile (which learning profile, color or large); "block"; 1 = colorful, 2 = large
%  3 - trial number
%  4 - resets for learning profile
%  5 - answer other (participant rating)
%  6 - outcome other (feedback)
%  7 - answer self (color rating for each object)
%  8 - answer self (large rating for each object)
%  9 - mean (average rating for color for each item)
% 10 - mean (average rating for large for each item)
% 11 - median (median rating for color for each item)
% 12 - median (median rating for large for each item)
% 13 - item list (item number)
% 14 - primary category
% 15 - secondary category 
% 16 - reset in subcat

%% load data

load ('./02_Main_experiment_computational_modeling/Data/nsDataSort_missing_trials_removed.mat');

%% fix reset

numSub=size(nsDataSort.data_all_scat,2);
nsDataSort.data_all={};
nsDataSort.data_all_cat={};

%% here add new item number that matches with size of sweep corr and add reset vector
% trial sort
for i_sub=1:numSub
n=size(nsDataSort.data_all_scat{1,i_sub},1) ;
nsDataSort.data_all{1,i_sub}= sortrows(nsDataSort.data_all_scat{1,i_sub},[6,4]); 
end

% reset var
 for i_sub=1:numSub
 resetVal = [abs(diff(nsDataSort.data_all{1,i_sub}{:,4} ))];
 resetVal(resetVal<10)=0;
 resetVal(resetVal>10)=1;
 resetVal =[1; resetVal];
 nsDataSort.data_all{1,i_sub}{:,7}=resetVal; 
 end

% cat sort
for i_sub=1:numSub
n=size(nsDataSort.data_all_scat{1,i_sub},1) ;
nsDataSort.data_all_cat{1,i_sub}= sortrows(nsDataSort.data_all_scat{1,i_sub},[6,15,4]); 
end

 for i_sub=1:numSub
 resetVal1 = [1; diff(nsDataSort.data_all_cat{1,i_sub}{:,15} )];
 resetVal1(resetVal1~=0)=1;
 nsDataSort.data_all_cat{1,i_sub}{:,7}=resetVal1; 
 end

% sort subcat
for i_sub=1:numSub
nsDataSort.data_all_scat{1,i_sub} = sortrows(nsDataSort.data_all_scat{1,i_sub},[6,16,4]);
end

% subcat sort
for i_sub=1:numSub
resetVal2 = [1; diff(nsDataSort.data_all_scat{1,i_sub}{:,16} )];
resetVal2(resetVal2 ~= 0) = 1;
nsDataSort.data_all_scat{1,i_sub}{:,7}=resetVal2; 
end

%% save data for modeling

save('/Volumes/research/Rosenblau/ANDI_data/Nonsocial_learning_task/modeling/data/nsDataSortModel.mat','nsDataSort' );

