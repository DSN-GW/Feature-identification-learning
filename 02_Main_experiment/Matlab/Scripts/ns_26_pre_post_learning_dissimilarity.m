%% Pre/post learning dissimilarity
% Set directories
matpath = matlab.desktop.editor.getActiveFilename;
idx     = strfind(matpath,'\');
dr.top  = matpath(1:idx(end-3)-1);
dr.exp  = '/02_Main_experiment/Matlab';
dr.data = '/Data/Raw_data';
dr.rpro = [dr.top,'/02_Main_experiment/R/Data'];
addpath(genpath(dr.top));
cd([dr.top,dr.exp,dr.data]);
clear idx

% Load data
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
load('item_ind_all.mat') % Stimulus info
exp1 = load('allrdms_exp1.mat');
exp2 = load('avg_rdms_exp2.mat');
load('arranged_names_exp2.mat');
exp2.noms = noms;
clear noms

%% From exp1, grab only items used in exp2
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
% Average exp1 rdms
exp1.avg_rdm = squareform(nanmean(exp1.allrdms,2));

% Get only activities and food
idx = item_ind.cat1 == 1 | item_ind.cat1 == 3;
exp1.rdm_no_fashion = exp1.avg_rdm(idx,idx);
exp1.noms = item_ind.filename(item_ind.cat1 == 1 | item_ind.cat1 == 3);

% Find items in exp1 that were used in exp2
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
% Color
idx = ismember(string(exp1.noms),string(exp2.noms.color));
exp1.color = exp1.rdm_no_fashion(idx,idx);

% Visualize
figpos = [10 10 800 700];
figure('Position',figpos)
imagesc(exp1.color); colorbar
xticks(1:42)
yticks(1:42)
xticklabels(exp1.noms(idx))
yticklabels(exp1.noms(idx))
xtickangle(90)
title('experiment 1 dissimilarities (colorful profile)')
set(gca,'FontSize',12)

% Large
idx = ismember(string(exp1.noms),string(exp2.noms.large));
exp1.large = exp1.rdm_no_fashion(idx,idx);

% Visualize
figpos = [10 10 800 700];
figure('Position',figpos)
imagesc(exp1.large); colorbar
xticks(1:42)
yticks(1:42)
xticklabels(exp1.noms(idx))
yticklabels(exp1.noms(idx))
xtickangle(90)
title('experiment 1 dissimilarities (large profile)')
set(gca,'FontSize',12)

%% Exp2 rdms
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
% Visualize color
exp2.avg.color = squareform(exp2.avg.color);
figpos = [10 10 800 700];
figure('Position',figpos)
imagesc(exp2.avg.color); colorbar
xticks(1:42)
yticks(1:42)
xticklabels(string(exp2.noms.color))
yticklabels(string(exp2.noms.color))
xtickangle(90)
title('experiment 2 dissimilarities (colorful profile)')
set(gca,'FontSize',12)

% Visualize large
exp2.avg.large = squareform(exp2.avg.large);
figpos = [10 10 800 700];
figure('Position',figpos)
imagesc(exp2.avg.large); colorbar
xticks(1:42)
yticks(1:42)
xticklabels(string(exp2.noms.large))
yticklabels(string(exp2.noms.large))
xtickangle(90)
title('experiment 2 dissimilarities (large profile)')
set(gca,'FontSize',12)

%% Differences
% Color
diffs.color = exp2.avg.color - exp1.color;
figpos = [10 10 800 700];
figure('Position',figpos)
imagesc(diffs.color); colorbar
xticks(1:42)
yticks(1:42)
xticklabels(string(exp2.noms.color))
yticklabels(string(exp2.noms.color))
xtickangle(90)
title('dissimilarity differences (colorful profile)')
set(gca,'FontSize',12)

% Large
diffs.large = exp2.avg.large - exp1.large;
figpos = [10 10 800 700];
figure('Position',figpos)
imagesc(diffs.large); colorbar
xticks(1:42)
yticks(1:42)
xticklabels(string(exp2.noms.large))
yticklabels(string(exp2.noms.large))
xtickangle(90)
title('dissimilarity differences (large profile)')
set(gca,'FontSize',12)

%% Visualizing both exps and diffs in one figure
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
cup1 = [-4,1];
cup2 = [-.5,1];
figpos = [10 10 1400 700];
figure('Position',figpos)
subplot(2,3,1)
h = imagesc(zscore(exp1.color)); colorbar
caxis(cup1)
ylabel('colorful profile')
xticks([])
yticks([])
title('experiment 1 dissimilarities')
set(gca,'FontSize',12)
subplot(2,3,4)
imagesc(zscore(exp1.large)); colorbar
caxis(cup1)
xticks([])
yticks([])
ylabel('large profile')
set(gca,'FontSize',12)
subplot(2,3,2)
imagesc(zscore(exp2.avg.color)); colorbar
caxis(cup1)
title('experiment 2 dissimilarities')
xticks([])
yticks([])
set(gca,'FontSize',12)
subplot(2,3,5)
imagesc(zscore(exp2.avg.large)); colorbar
caxis(cup1)
xticks([])
yticks([])
set(gca,'FontSize',12)
df1 = zscore(exp1.color) - zscore(exp2.avg.color);
df2 = zscore(exp1.large) - zscore(exp2.avg.large);
subplot(2,3,3)
imagesc(df1); colorbar
caxis(cup2)
title('differences')
xticks([])
yticks([])
set(gca,'FontSize',12)
subplot(2,3,6)
imagesc(df2); colorbar
caxis(cup2)
xticks([])
yticks([])
set(gca,'FontSize',12)

%% add difference matrices into regression df
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
df = readtable('df_regression_exp2.csv');
df.Var7 = [squareform(diffs.color)';squareform(diffs.large)'];
writetable(df,'df_regression_exp2.csv');

%% save structure for permutation test
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
df_perm = [squareform(exp1.color)';squareform(exp1.large)';
           squareform(exp2.avg.color)';squareform(exp2.avg.large)'];
df_perm = [df_perm,[repmat(0,size(df_perm,1)/2,1);repmat(1,size(df_perm,1)/2,1)]];
condition = [repmat(0,size(df_perm,1)/2/2,1);repmat(1,size(df_perm,1)/2/2,1)];
df_perm = [df_perm,[condition;condition]];
writematrix(df_perm,'df_permutation_testing.csv');

%% Demean/normalize within each sample
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
% subtract mean from all dissimilarities for each experiment
% then take absolute difference
% then do exp2 - exp1 dissimilarities
% both should be positive (because of absolute)
% sign is meaningful such that negative meanse exp2 is smaller than exp1
format longg
mean_exp1.clr = mean(squareform(exp1.color));
mean_exp1.lrg = mean(squareform(exp1.large));
mean_exp2.clr = mean(squareform(exp2.avg.color));
mean_exp2.lrg = mean(squareform(exp2.avg.large)); % very close, but diff.

% Demean
exp1_demeaned.clr = abs(squareform(exp1.color)' - mean_exp1.clr);
exp1_demeaned.lrg = abs(squareform(exp1.large)' - mean_exp1.lrg);
exp2_demeaned.clr = abs(squareform(exp2.avg.color)' - mean_exp2.clr);
exp2_demeaned.lrg = abs(squareform(exp2.avg.large)' - mean_exp2.lrg);

% Find limits
lims = zeros(6,2,1);
lims(1,:) = [min(exp1_demeaned.clr),max(exp1_demeaned.clr)];
lims(2,:) = [min(exp1_demeaned.lrg),max(exp1_demeaned.lrg)];
lims(3,:) = [min(exp2_demeaned.clr),max(exp2_demeaned.clr)];
lims(4,:) = [min(exp2_demeaned.lrg),max(exp2_demeaned.lrg)];
lims(5,:) = [min(exp2_demeaned.clr - exp1_demeaned.clr),max(exp2_demeaned.clr - exp1_demeaned.clr)];
lims(6,:) = [min(exp2_demeaned.lrg - exp1_demeaned.lrg),max(exp2_demeaned.lrg - exp1_demeaned.lrg)];

% Min and max
min(lims,[],"all")
max(lims,[],"all")

% Visualize
clim = [min(lims,[],"all"),.0085];
tiledlayout(3,2,'TileSpacing','compact')
nexttile(1)
imagesc(squareform(exp1_demeaned.clr)); colorbar
caxis(clim)
ylabel('exp1 dissimilarities (demeaned)')
title('color profile','FontWeight','light')
nexttile(2)
imagesc(squareform(exp1_demeaned.lrg)); colorbar
caxis(clim)
title('large profile','FontWeight','light')
nexttile(3)
imagesc(squareform(exp2_demeaned.clr)); colorbar
caxis(clim)
ylabel('exp2 dissimilarities (demeaned)')
nexttile(4)
imagesc(squareform(exp2_demeaned.lrg)); colorbar
caxis(clim)
nexttile(5)
imagesc(squareform(exp2_demeaned.clr - exp1_demeaned.clr)); colorbar
caxis(clim)
ylabel('exp2 - exp1 dissimilarities (demeaned)')
nexttile(6)
imagesc(squareform(exp2_demeaned.lrg - exp1_demeaned.lrg)); colorbar
caxis(clim)
diffs.demeaned.clr = exp2_demeaned.clr - exp1_demeaned.clr;
diffs.demeaned.lrg = exp2_demeaned.lrg - exp1_demeaned.lrg;

% Add vars to df
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
df = readtable('df_regression_exp2.csv');
df.Var8 = [diffs.demeaned.clr;diffs.demeaned.lrg];
df.Var9 = [exp2_demeaned.clr;exp2_demeaned.lrg];
df.Var10 = [exp1_demeaned.clr;exp1_demeaned.lrg];
writetable(df,'df_regression_exp2.csv');




