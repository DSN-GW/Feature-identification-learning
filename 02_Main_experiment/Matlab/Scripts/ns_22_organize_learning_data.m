%% Organize learning data
clear
close all

% Set directories
matpath = matlab.desktop.editor.getActiveFilename;
idx     = strfind(matpath,'\');
dr.top  = matpath(1:idx(end-3)-1);
dr.exp  = '/02_Main_experiment/Matlab';
dr.data = '/Data/Raw_data';
addpath(genpath(dr.top));
cd([dr.top,dr.exp,dr.data]);
clear idx

% Load data
% note that pre-experiment: exp1 and main experiment: exp2 in naming
load('item_ind_all.mat')                    % Stimulus info
exp1 = load('full_data_struct_exp1.mat');   % Exp1 data
exp2 = load('data_structure_full_n92.mat'); % Exp2 data
ref_data = load('nsDataSotModel.mat');      % From GR

%% Structure reference
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

%% First extract rating and compute means/medians from experiment 1
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
% Note: we will not use expensive ratings in our analyses, so we skip them
% Note: we also only used a subset of stimuli (84/120)
% Note: ratings included new stimuli (not used in the experiment), these
% were replacements for the fashion items. because of this, we need to be
% mindful of filenames across structures.
redcap.color = [exp1.data.redcap.data.subid,exp1.data.redcap.data(:,31:3:390)];
redcap.large = [exp1.data.redcap.data.subid,exp1.data.redcap.data(:,33:3:390)];
redcap.color.Properties.VariableNames{1} = 'subid';
redcap.large.Properties.VariableNames{1} = 'subid';

% Compute means/medians, note that items (columns) are sorted alphabetically
means.color = varfun(@mean, redcap.color(:,2:end), 'InputVariables', @isnumeric);
means.large = varfun(@mean, redcap.large(:,2:end), 'InputVariables', @isnumeric);
medians.color = varfun(@median, redcap.color(:,2:end), 'InputVariables', @isnumeric);
medians.large = varfun(@median, redcap.large(:,2:end), 'InputVariables', @isnumeric);

% Create matrix with all relevant info
exp1_info =zeros(120,4);
exp1_info(:,1) = (table2array(varfun(@mean, redcap.color(:,2:end), 'InputVariables', @isnumeric)))';
exp1_info(:,2) = (table2array(varfun(@mean, redcap.large(:,2:end), 'InputVariables', @isnumeric)))';
exp1_info(:,3) = (table2array(varfun(@median, redcap.color(:,2:end), 'InputVariables', @isnumeric)))';
exp1_info(:,4) = (table2array(varfun(@median, redcap.large(:,2:end), 'InputVariables', @isnumeric)))';

% Now remove prefix from item names
for i_col = 1:size(means.color,2)
    means.color.Properties.VariableNames{i_col} = ...
        means.color.Properties.VariableNames{i_col}(10:end);
    means.large.Properties.VariableNames{i_col} = ...
        means.large.Properties.VariableNames{i_col}(10:end);
    medians.color.Properties.VariableNames{i_col} = ...
        medians.color.Properties.VariableNames{i_col}(12:end);
    medians.large.Properties.VariableNames{i_col} = ...
        medians.large.Properties.VariableNames{i_col}(12:end);
end

% Drop prefixes and transpose raw rating data for referencing later
for i_tem = 2:size(redcap.color,2)
    redcap.color.Properties.VariableNames{i_tem} = ...
        redcap.color.Properties.VariableNames{i_tem}(5:end);
    redcap.large.Properties.VariableNames{i_tem} = ...
        redcap.large.Properties.VariableNames{i_tem}(5:end);
end

% Sort item ind alphabetically
vars_redcap = redcap.color.Properties.VariableNames(2:end);
[ind1,ind2]=sort(vars_redcap);
exp1_info(:,5)=ind2';
ind1=string(ind1);%
item_ind = sortrows(item_ind,'filename');
sum(strcmp(item_ind.filename, (ind1)'))==120; % If 1 all good
exp1_info=sortrows(exp1_info,5); % sorted
item_ind1=array2table(exp1_info(:,1:4));
item_ind1.Properties.VariableNames={'mean_color','mean_large','med_color','med_large'};

% New matrix
item_ind_new=[item_ind item_ind1];

%% Now we can start building the final structure
% First, extract the learning data and initialize final structure
psiturk = exp2.data.psiturk;
nsDataSort.data_all_scat = cell(1,size(psiturk.subids,1));

% We also need raw rating data (i.e., self ratings) from experiment 2
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
redcap.color = [exp2.data.redcap.data.subid,exp2.data.redcap.data(:,36:3:395)];
redcap.large = [exp2.data.redcap.data.subid,exp2.data.redcap.data(:,38:3:395)];
redcap.color.Properties.VariableNames{1} = 'subid';
redcap.large.Properties.VariableNames{1} = 'subid';

%% Now, extract participant level data and place into structure
for i_sub = 1:size(psiturk.subids,1)
    % Participant data
    subid = psiturk.subids{i_sub};
    d = psiturk.data(strcmp(psiturk.data.subid,subid),:);

    % Because our word-concept pairings change, the pseudowords are not
    % useful for us. However, block is, with 1 = colorful and 2 = large
    d = d(:,[1,3:end]); % Drop unneeded 'profile' column
    d.Properties.VariableNames{6} = 'profile'; % Rename block

    % Add profile reset variable and placeholder variables
    d.reset = [1;abs(diff(d.profile))]; % Profile reset column
    d.self_color = nan(size(d,1),1);
    d.self_large = nan(size(d,1),1);
    d.mean_color = nan(size(d,1),1);
    d.mean_large = nan(size(d,1),1);
    d.median_color = nan(size(d,1),1);
    d.median_large = nan(size(d,1),1);
    d.itemnr_new = nan(size(d,1),1);
    d.cat1 = nan(size(d,1),1);
    d.cat2 = nan(size(d,1),1);
    d.reset_subcat = nan(size(d,1),1);
    d.PEsigned = nan(size(d,1),1);
    d.PEabsolute = nan(size(d,1),1);

    % Add back to overall structure
    nsDataSort.data_all_scat{i_sub} = d;
end

% Sort alphabetically
for i_sub = 1:size(psiturk.subids,1)
    nsDataSort.data_all_scat{i_sub} = sortrows(nsDataSort.data_all_scat{i_sub},'filename');
end

%% Here extract the item sorting from redcap
redcap_self_vars = [redcap.color.Properties.VariableNames;redcap.large.Properties.VariableNames];
redcap_self_vars(:,1) = [];

% Now remove prefix from item names
for i_col = 1:size(redcap_self_vars,2)
    for i_row = 1:2
        redcap_self_vars{i_row,i_col} = ...
            redcap_self_vars{i_row,i_col}(5:end);
    end
end
sum(strcmp(redcap_self_vars(1,:),redcap_self_vars(2,:)))==120; % 1 we are good

% Double check if items have the same sorting across learning task and
% redcap self; there are only 84 items so identify them
learning_items = nsDataSort.data_all_scat{i_sub}.filename;
learning_items = erase(learning_items,'.png');
learning_items_str = string(learning_items);
for i_red = 1:length(redcap_self_vars)
    for i_ler=1:length(learning_items)
        cmp(i_ler,i_red)=  strcmp(learning_items(i_ler),redcap_self_vars(1,i_red));
        cmp2(i_ler,i_red)= strcmp(learning_items_str(i_ler),item_ind_new.filename(i_red,1));
    end
end
select = sum(cmp,1);
select2 = sum(cmp2,1);
vars_select=[1 select];

% Take from redcap only vars that are in common
redcap.color_new=redcap.color(:,vars_select==1);
redcap.large_new=redcap.large(:,vars_select==1);

% Take from exp 1 vars in common
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
item_ind_new=item_ind_new(find(select2==1)',:);

%% Now add everything into the preexisting table
for i_sub = 1:size(psiturk.subids,1)
    row = find(strcmp(psiturk.subids{i_sub},redcap.color_new.subid));
    try
        nsDataSort.data_all_scat{i_sub}.self_color = (redcap.color_new{row,2:end})';
        nsDataSort.data_all_scat{i_sub}.self_large = (redcap.large_new{row,2:end})';
        nsDataSort.data_all_scat{i_sub}.mean_color = item_ind_new.mean_color;
        nsDataSort.data_all_scat{i_sub}.mean_large = item_ind_new.mean_large;
        nsDataSort.data_all_scat{i_sub}.median_color = item_ind_new.med_color;
        nsDataSort.data_all_scat{i_sub}.median_large = item_ind_new.med_large;
        nsDataSort.data_all_scat{i_sub}.itemnr_new = item_ind_new.itemnr_new;
        nsDataSort.data_all_scat{i_sub}.cat1 = item_ind_new.cat1;
        nsDataSort.data_all_scat{i_sub}.cat2 = item_ind_new.cat2;
    catch
        fprintf('Inconsistent data in iteration %s, skipped.\n', i);
        clear row
    end
end

%% Adding PE variables and saving
for i_sub = 1:size(nsDataSort.data_all_scat,2)
    % compute signed PE
    nsDataSort.data_all_scat{i_sub}.PEsigned = nsDataSort.data_all_scat{i_sub}.feedback - ...
        nsDataSort.data_all_scat{i_sub}.rating;
    % compute absolute PE
    nsDataSort.data_all_scat{i_sub}.PEabsolute = abs(nsDataSort.data_all_scat{i_sub}.PEsigned);
    % add resets for subcategory transitions
    nsDataSort.data_all_scat{i_sub}.reset_subcat = [1;diff(nsDataSort.data_all_scat{i_sub}.cat2)>0];
end
cd([dr.top,dr.exp,'/Data'])
save('nsDataSort_before_exclusions.mat','nsDataSort')

%% Next, compute proportions of missing data to make exclusion decisions
propMissing = cell(3,size(nsDataSort.data_all_scat,2)); % Initialize
for i_sub = 1:size(nsDataSort.data_all_scat,2)
    % Get participant data
    d = nsDataSort.data_all_scat{i_sub};
    subid = d.subid{1};

    % Nonresponse trials are coded as zero
    numMissing = sum(d.rating==0);

    % Save id and proportion
    propMissing{1,i_sub} = subid;
    propMissing{2,i_sub} = numMissing;
    propMissing{3,i_sub} = numMissing/size(d,1);

    % Mark nonresponse trials as NaN for removal later
    d.rating(d.rating==0) = NaN;
    nsDataSort.data_all_scat{i_sub} = d;
end

% Visualize distribution
figure
subplot(1,2,1)
histogram(cell2mat(propMissing(2,:)))
xlabel('num missing trials')
ylabel('frequency')
subplot(1,2,2)
histogram(cell2mat(propMissing(3,:)))
xlabel('proportion missing trials')

% Use mean + 2SD as cutoff for participant exclusions
m = mean(cell2mat(propMissing(2,:)));
s = std(cell2mat(propMissing(2,:)));
c = m+(2*s);
e = cell2mat(propMissing(2,:))>c;  %index of exclusions
sum(cell2mat(propMissing(2,:))>c); %number of exclusions
ids_to_exclude = propMissing(1,e)';

%% exclude participants with excessive missing data & save structs
nsDataSort.data_all_scat = nsDataSort.data_all_scat(1,~e); % Participants

% Check for participants with too many trials
sessions = []; % To store session lengths
for i_sub = 1:size(nsDataSort.data_all_scat,2)
    sessions(i_sub) = size(nsDataSort.data_all_scat{i_sub},1);
end
histogram(sessions)
idx = sessions > 84;
checking = nsDataSort.data_all_scat(idx); % Which participants?

%% Exclusions
% Note: investigating which trials can be kept and which can be excluded
% for those participants that have excess trials. There are 5 participants 
% with more than 84 trials. Extra trials range from just one extra trial to
% another entire learning block. For those that clearly completed an extra 
% thrid block, and for those who had a false start (i.e., one or so extra 
% trials at the start, we will just drop the extra trials. For those that 
% completed several trials and then restarted, their data are likely 
% contaminated and cannot be used, so we will exclude them entirely.

% Participant {'613bb7f2a06a3b358fe33f46'} completed 126 trials with
% profiles ordered 2, 1, 1. so, we will just drop the third block of trials
% indices = find(idx); % Find participant indices
% nsDataSort.data_all_scat{indices(1)} = nsDataSort.data_all_scat{indices(1)}(1:84,:);

% Participant {'5e95116fcd40fa88e616d9d9'} completed 23 trials in profile 1
% before starting over and completing full blocks of profile 2, then 1. So,
% we must decide whether the first 23 trials of profile 1 contaminate the
% later trials of profile 1 (maybe even profile 2). Best to exclude
% altogether. We'll grab the index here, then exclude after checking others
% excl = indices(2);

% Participant {'61085282a230691c4078e3c2'} completed 96 trials with one
% complete block of profile 1 first, then 12 trials of profile 2 before
% starting over with a full block of profile 2. The twelve extra trials
% immediately preceeding the full block are likely contaminanted; exculde
% excl = [excl,indices(3)];

% Participant {'61098efe1c7dcf4133162908'} completed 85 trials (1 more than
% the correct total). This extra trial occurred at the start of the session.
% While the item seen on the first trial is repeated on a subsequent trial
% (23), a single exposure should not impact the results; drop extra trial
% nsDataSort.data_all_scat{indices(4)} = nsDataSort.data_all_scat{indices(4)}(2:end,:);

% Participant {'6155bb4f5250fc58f16eca99'} completed 127 trials (i.e., 3
% full blocks + 1 extra trial). One extra trial occurred at the start, like
% the last participant. They then went to complete 2 full blocks of profile
% 2 and one full block of profile 1. Because the blocks are intact, we can
% remove the second block of profile 2 as it is certainly contaminated by
% the first block. We can retain the full block of profile 1 since it should
% not be contaminated by the extra exposures to profile 2.
% nsDataSort.data_all_scat{indices(5)} = nsDataSort.data_all_scat{indices(5)}([2:43,86:end],:);

% Save all excluded ids for later reference
% ids_to_exclude = [ids_to_exclude;{'5e95116fcd40fa88e616d9d9'};{'61085282a230691c4078e3c2'}];


%% Now make final exclusions
% idx = 1:size(nsDataSort.data_all_scat,2);
% idx = ismember(idx,excl);
% nsDataSort.data_all_scat = nsDataSort.data_all_scat(~idx);
indices = find(idx); % Find participant indices
idx = 1:size(nsDataSort.data_all_scat,2);
idx = ismember(idx,indices);
nsDataSort.data_all_scat = nsDataSort.data_all_scat(~idx);

%% Gender
gender = table(exp2.data.redcap.data.subid,exp2.data.redcap.data.gender);
for i = 1:size(nsDataSort.data_all_scat,2)
    for j = 1:size(gender,1)
        if strcmp(nsDataSort.data_all_scat{i}.subid{1},gender.Var1{j})
            nsDataSort.data_all_scat{i}.gender = repmat(gender.Var2(j),size(nsDataSort.data_all_scat{i},1),1);
        end
    end
end
cd([dr.top,dr.exp,'/Data'])
save('nsDataSort_after_exclusions.mat','nsDataSort')

%% And loop through structure to remove missing trials for later modeling
for i_sub = 1:size(nsDataSort.data_all_scat,2)
    nsDataSort.data_all_scat{i_sub} = rmmissing(nsDataSort.data_all_scat{i_sub});
end
cd([dr.top,dr.exp,'/Data'])
save('nsDataSort_missing_trials_removed.mat','nsDataSort')
