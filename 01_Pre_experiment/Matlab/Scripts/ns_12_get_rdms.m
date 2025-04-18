%% Generates RDM structures
clear
close all

% Set directories
matpath = matlab.desktop.editor.getActiveFilename;
idx     = strfind(matpath,'\');
dr.top  = matpath(1:idx(end-2)-1);
dr.exp  = '/01_Pre_experiment';
dr.data = '/Data/Raw_data';
addpath(genpath(dr.top));
cd([dr.top,dr.exp,dr.data]);

% Read in stimulus info and full data structure
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
load('full_data_struct_exp1.mat'); % data
load('item_ind_all.mat'); % stimulus info


%% Correcting filename discrepancy
%  Due to filename mismatch between recorded meadows data and item_ind,
%  necessary to correct filenames after extracting from meadows. This is
%  done in the main loop before arranging, since the arrange_rdm function
%  extracts item numbers by matching filenames. Corrections are listed here
%  where column 1 is the extracted name and column 2 is the corrected name
file_corrections = {'accordion_4'   ,      'accordion_2'     ;
                    'apple_02'      ,      'apple_2'         ;
                    'apple_04'      ,      'apple_4'         ;
                    'berries_5'     ,      'berries_05'      ;
                    'bowling_3'     ,      'bowling_03'      ;
                    'chocolate2'    ,      'chocolate_2'     ;
                    'cologne_01'    ,      'cologne_1'       ;
                    'cologne_02'    ,      'cologne_2'       ;
                    'green_beans01' ,      'green_beans_1'   ;
                    'greens_05'     ,      'greens_07'       ;
                    'makeup_04'     ,      'makeup_4'};
                
%% Main processing loop
% Initialize vars
newdata = data.meadows.data;
allrdms = [];
all_with_cats = [];
ind = 0;
rdm_mats={};
ntrials = [];

% Loop through subjects
for i_sub = 1:size(newdata,1)
    % Under rdm vector of item-pair dissimilarity// 120 items
    % n*(n-1)/2 = 7140 pairs
    % needed for correctly indexing into meadows structure
    if i_sub < 77 || i_sub > 85
        ntask = 2; % One is instruction and one is data
    else
        ntask = 3; % Here first 2 structures are instruction and third is data
    end
    
    % Extract rdm and stimulus names
    rdm   = newdata{i_sub}.tasks{ntask}.rdm;
    names = struct2cell(newdata{i_sub}.tasks{ntask}.stimuli);
    names = names(2,:)';
    ntrials(i_sub,1) = size(newdata{i_sub}.tasks{ntask}.trials,1);
    
    % Correct filenames before arranging
    for i_cor = 1:size(file_corrections,1)
        for i_tem = 1:size(names,1)
            if strcmp(file_corrections{i_cor,1},names{i_tem})
                names{i_tem} = file_corrections{i_cor,2};
            end
        end
    end
    
    % Arrange rdm by item number and category info
    [arranged_rdm, arranged_names] = arrange_rdm(rdm,names,item_ind);
    rdm_mats{i_sub,1}.data = arranged_rdm; %120x120 mat
    rdm_mats{i_sub,1}.sort = arranged_names; %order file for each participant
    
    % Aggregate rdms 
    allrdms = [allrdms, squareform(arranged_rdm)'];
    clear arranged_rdm arranged_names
end    

%% save data structures 
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
save(fullfile([dr.top,dr.exp,'/data/', 'allrdms_exp1.mat']),'allrdms');
save(fullfile([dr.top,dr.exp,'/data/', 'rdm_sqmat.mat']),'rdm_mats'); % saves sort and matrix 

%% Visualize average rdm
% Average over subject RDMs
avgrdm = nanmean(allrdms,2);

% For axis labels below
for i_chr = 1:size(item_ind,1)
arranged_names{i_chr} = char(item_ind.filename{i_chr});
end

% Visualize
figpos = [10 10 1100 1000];
figure('Position',figpos)
h = imagesc(squareform(avgrdm./nanmax(avgrdm,[],'all'))); colorbar
xticks([])
% xticklabels(arranged_names)
% xtickangle(90)
yticks([])
% yticklabels(arranged_names)
set(gca,'FontSize',20)
title('Average RDM','FontSize',20)

%% Check average mat from indiv structure
matrdm(120,120,245)=NaN;
for i_sub=1:size(newdata,1)
matrdm(:,:,i_sub) = rdm_mats{i_sub,1}.data;
end
figure
imagesc (nanmean(matrdm,3));
    
%% Visualize correlation matrix
% Overall from study 1
figpos = [10 10 1100 1000];
figure('Position',figpos)
rdm = squareform(avgrdm);
imagesc(corr(rdm)); colorbar
xticks(1:120)
xticklabels(arranged_names)
xtickangle(90)
yticks(1:120)
yticklabels(arranged_names)
set(gca,'FontSize',8)
title('average rdm (correlation matrix)','FontSize',20)
sweep_corr_s1_full = corr(rdm);

% Save data
save('sweep_corr_s1_full.mat','sweep_corr_s1_full')

%% Construct category matrix for regression modeling
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
cat1_vector = item_ind.cat1;
cat2_vector = item_ind.cat2;
item_means = mean(squareform(avgrdm),2);
df = [item_means,cat1_vector,cat2_vector];
csvwrite('df_regression_exp1.csv',df)










