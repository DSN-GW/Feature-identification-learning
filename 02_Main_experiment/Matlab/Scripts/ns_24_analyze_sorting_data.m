%% Analyze sorting data
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
load('item_ind_all.mat') % Stimulus info
load('data_structure_full_n92.mat')

%% First, exclude ids identified in learning analysis
% Subjects had too much missing data (n=5, m+2sd), or too many trials (n=5)
ids.exclude = [{'5caa8efdf9a1bd0001ada0de';'60fe3799642da04fc85c5dbc';...
    '6110a3cea048e6a387686a00';'613927703e23f814407dc5c7';...
    '61481042c3f4bbe5169303e8';'5e95116fcd40fa88e616d9d9';...
    '61085282a230691c4078e3c2';'613bb7f2a06a3b358fe33f46';
    '61098efe1c7dcf4133162908';'6155bb4f5250fc58f16eca99'}];
% Find exclusion indices
for i_sub = 1:size(data.meadows.color.ids,1)
    idx.color(i_sub) = ismember(data.meadows.color.ids{i_sub,2},ids.exclude);
    idx.large(i_sub) = ismember(data.meadows.large.ids{i_sub,2},ids.exclude);
end
% Remove
data.meadows.color.ids = data.meadows.color.ids(~idx.color,:);
data.meadows.large.ids = data.meadows.large.ids(~idx.large,:);

% Save data
cd([dr.top,dr.exp,'/Data'])
save('data_structure_full_n82.mat','data')


%% Extract rdms
% Initialize vars
rdms.color = [];
rdms.large = [];
ntrials.color = [];
ntrials.large = [];
nanper.color = [];
nanper.large = [];
gnd = {};

% Loop through colorful and large conditions
for i_blk = 1:size(fieldnames(data.meadows),1)

    % Initialize
    ntri = [];
    nancnt = [];

    % Set condition variables
    if i_blk == 1
        condition = data.meadows.color;
        condtag = 'colorful';
    elseif i_blk == 2
        condition = data.meadows.large;
        condtag = 'large';
    end

    % Gather rdms
    % Note that pre-experiment: exp1 and main experiment: exp2 in naming
    sub_rdms = [];
    for i_sub = 1:size(condition.ids,1)
        gnd{i_sub,1} = condition.ids{i_sub,2};
        rdm = condition.ids{i_sub,1}.tasks{2,1}.rdm';
        filenames = struct2cell(condition.ids{i_sub,1}.tasks{2,1}.stimuli);
        filenames = filenames(2,:)';
        [arranged_rdm, arranged_names, with_cats] = arrange_rdm_exp2(rdm,filenames,item_ind);
        arranged_rdm = squareform(arranged_rdm);
        sub_rdms = [sub_rdms;arranged_rdm];
        ntri(i_sub,1) = size(condition.ids{i_sub,1}.tasks{2,1}.trials,1);
        nancnt(i_sub,1) = sum(isnan(rdm))/length(rdm);
    end
    clear rdm

    % Assign to correct structure
    if i_blk == 1
        rdms.color = sub_rdms;
        noms.color = arranged_names;
        ntrials.color = ntri;
        nanper.color = nancnt;
    elseif i_blk == 2
        rdms.large = sub_rdms;
        noms.large = arranged_names;
        ntrials.large = ntri;
        nanper.large = nancnt;
    end
end


%% Average rdms for each condition, then visualize
for i_tem = 1:size(noms.color,1)
    noms.color(i_tem) = strrep(noms.color{i_tem},'_','');
    noms.large(i_tem) = strrep(noms.large{i_tem},'_','');
end
ind = zeros(41,1);
ind(1:2:41) = 1;
ind = [ind(2:end);0;1];

% Compute means
avg.color = mean(rdms.color);
avg.large = mean(rdms.large);

% Visualize
cup = [0,0.04];
figpos = [10 10 1800 500];
figure('Position',figpos)
subplot(1,3,1)
imagesc(squareform(avg.color)); colorbar
ax = gca;
ax.CLim = cup;
yticks(1:2:42)
yticklabels(string(noms.color(ind==1)))
xticks(1:2:42)
xticklabels(string(noms.color(ind==1)))
xtickangle(90)
title('colorful')
set(gca,'FontSize',15)
subplot(1,3,2)
imagesc(squareform(avg.large)); colorbar
ax = gca;
ax.CLim = cup;
yticks(1:2:42)
yticklabels(string(noms.large(ind==1)))
xticks(1:2:42)
xticklabels(string(noms.large(ind==1)))
xtickangle(90)
title('large')
set(gca,'FontSize',15)
% Just to check that the values are different
subplot(1,3,3)
avg.diff = avg.color - avg.large;
imagesc(squareform(avg.diff)); colorbar
title('difference')
xlabel('item number')
ylabel('item number')
set(gca,'FontSize',15)

% Save structures
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
cd([dr.top,dr.exp,'/Data'])
save('avg_rdms_exp2.mat','avg')
save('arranged_names_exp2.mat','noms')

%% Checking for gender differences in dissimilarities
for i_sub = 1:length(gnd)
    for i_row = 1:length(data.redcap.data.subid)
        if strcmp(gnd{i_sub,1},data.redcap.data.subid{i_row})
            gnd{i_sub,2} = data.redcap.data.gender(i_row);
        end
    end
end
gnd_col = [];
for i_sub = 1:size(gnd,1)
    gnd_col = [gnd_col;repmat(gnd{i_sub,2},size(rdms.color,2),1)];
end
tmp_col = rdms.color';
tmp_lrg = rdms.large';
dissims = [tmp_col(:);tmp_lrg(:)];
df_gnd = table(repmat(gnd_col,2,1),[repmat(1,length(gnd_col),1);repmat(2,length(gnd_col),1)],dissims);
df_gnd.Properties.VariableNames = {'gender','condition','dissim'}
cd([dr.top,dr.exp,'/Data'])
writetable(df_gnd,'df_gender_dissim.csv')

%% Arranging for bootstrapping exp2 dissimilarity estimates
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
%  GR's request in response to the reviewers
df = [[zeros(861,1); ones(861,1)], [rdms.color'; rdms.large']];
cd([dr.top,dr.exp,'/Data'])
writematrix(df,'df_rdms_exp2_for_bootstrapping.csv')

% Checking if bootstraps and means are similar
tmp1 = mean(df(:,2:end),2);
cd([dr.top,dr.exp,'/Data'])
writematrix(tmp1,'exp2_mean_before_bootstrapping.csv')