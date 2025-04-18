%% Prep data for regression analyses
clear
close all

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
% note that pre-experiment: exp1 and main experiment: exp2 in naming
load('item_ind_all.mat') % Stimulus info
load('avg_rdms_exp2.mat')
load('full_data_struct_exp1.mat')
noms_arranged = load('arranged_names_exp2.mat');

%% Extract ratings for each item and compute item rating averages
% Transpose tables, column 1 = items, 2:end = subjects
rating_color = rows2vars(data.redcap.data(:,31:3:390));
rating_large = rows2vars(data.redcap.data(:,33:3:390));
noms = rating_color{:,1}; %get item names
for i_nom = 1:size(noms,1)
    noms{i_nom} = noms{i_nom}(5:end); %drop prefix
end

% Gather ratings
ratings.color = rating_color{:,2:end};
ratings.large = rating_large{:,2:end};

% Compute means
m_color = mean(ratings.color,2);
m_large = mean(ratings.large,2);

% Add these values into item_ind
item_ind.color = zeros(120,1);
item_ind.large = zeros(120,1);
for i_tem = 1:size(item_ind,1)
    for i_rat = 1:size(m_color,1)
        if strcmp(item_ind.filename{i_tem},noms{i_rat})
            item_ind.color(i_tem) = m_color(i_rat);
            item_ind.large(i_tem) = m_large(i_rat);
        end
    end
end

%% Add rating info to rdms
load('data_structure_full_n82.mat')
clear noms
for i_blk = 1:size(fieldnames(data.meadows),1)
    sub_rdms = [];
    if i_blk == 1
        condition = data.meadows.color;
        condtag = 'colorful';
    elseif i_blk == 2
        condition = data.meadows.large;
        condtag = 'large';
    end
    
    % Only need one from each condition
    % Note that pre-experiment: exp1 and main experiment: exp2 in naming
    rdm = condition.ids{1,1}.tasks{2,1}.rdm';
    filenames = struct2cell(condition.ids{1,1}.tasks{2,1}.stimuli);
    filenames = filenames(2,:)';
    [arranged_rdm, arranged_names, with_cats_rats] = arrange_rdm_full_exp2(rdm,filenames,item_ind);
    
    % Set to correct structure
    if i_blk == 1
        color.with_cat_rats = with_cats_rats;
        color.arranged_names = arranged_names;
    elseif i_blk == 2
        large.with_cat_rats = with_cats_rats;
        large.arranged_names = arranged_names;
    end
end

%% Set up category/rating matrices
% Category
matrix_category = ones(42,42);
matrix_category(1:21,1:21) = 0;
matrix_category(22:end,22:end) = 0;
imagesc(matrix_category); colorbar

% Subcategory
c_sub = [];
l_sub = [];
for i_tem = 1:size(color.arranged_names,1)
    for i_row = 1:size(item_ind,1)
        if strcmp(color.arranged_names{i_tem},item_ind.filename{i_row})
            c_sub(i_tem,1) = item_ind.cat2(i_row);
        end
        if strcmp(large.arranged_names{i_tem},item_ind.filename{i_row})
            l_sub(i_tem,1) = item_ind.cat2(i_row);
        end
    end
end

% Generate matrices
% Color
mat_dim1 = repmat(c_sub,1,42);
mat_dim2 = repmat(c_sub',42,1);
matrix_subcategory_color = mat_dim1~=mat_dim2;
figure
imagesc(matrix_subcategory_color); colorbar
% Large
mat_dim1 = repmat(l_sub,1,42);
mat_dim2 = repmat(l_sub',42,1);
matrix_subcategory_large = mat_dim1~=mat_dim2;
figure
imagesc(matrix_subcategory_large); colorbar
% Disparity matrices
matrix_color = [];
matrix_large = [];
for i_row = 1:size(squareform(avg.color),1)
    for i_col = 1:size(squareform(avg.color),1)
        matrix_color(i_row,i_col) = color.with_cat_rats(i_row,46) - color.with_cat_rats(46,i_col);
        matrix_large(i_row,i_col) = large.with_cat_rats(i_row,47) - large.with_cat_rats(47,i_col);
    end
end

% Visualize rating disparity matrices
cup = [0,8];
figpos = [10 10 1800 700];
figure('Position',figpos)
subplot(1,2,1)
imagesc(abs(matrix_color)); colorbar
ax = gca;
ax.CLim = cup;
xticks(1:42)
xticklabels(string(color.arranged_names))
xtickangle(90)
yticks(1:42)
yticklabels(string(color.arranged_names))
title('absolute difference in color rating')
set(gca,'FontSize',12)
subplot(1,2,2)
imagesc(abs(matrix_large)); colorbar
ax = gca;
ac.CLim = cup;
xticks(1:42)
xticklabels(string(large.arranged_names))
xtickangle(90)
yticks(1:42)
yticklabels(string(large.arranged_names))
title('absolute difference in large rating')
set(gca,'FontSize',12)

% Set up df
df_regression = [avg.color',repmat(0,length(avg.color'),1),squareform(matrix_category)',squareform(matrix_subcategory_color)',squareform(matrix_color)',squareform(matrix_large)'; ...
                 avg.large',repmat(1,length(avg.color'),1),squareform(matrix_category)',squareform(matrix_subcategory_large)',squareform(matrix_color)',squareform(matrix_large)'];


%% Detour to pre/post learning dissimilarity script to add to df
ns_26_pre_post_learning_dissimilarity

%% Gather exp1 and exp2 ratings for permutation testing
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
% Experiment1 (pre-experiment)
noms = table2cell(rating_color(:,1));
for i_tem = 1:size(noms,1)
    noms{i_tem} = noms{i_tem}(5:end)
end
clr = table2array(rating_color(:,2:end));
clr = clr(:);
lrg = table2array(rating_large(:,2:end));
lrg = lrg(:);
df_ratings = table(repmat(noms,size(rating_color,2)-1,1),clr,lrg,zeros(size(clr,1),1));

% Experiment 2 (main experiment)
exp2_data = load('rating_data_redcap_exp2.mat');
clr2 = [];
lrg2 = [];
noms2 = exp2_data.rating_data{1,1}.Properties.RowNames;
for i_tem = 1:size(noms,1)
    noms2{i_tem} = noms2{i_tem}(5:end)
end
for i_sub = 1:size(exp2_data.rating_data,1)
    clr2 = [clr2;table2array(exp2_data.rating_data{i_sub,1})];
    lrg2 = [lrg2;table2array(exp2_data.rating_data{i_sub,3})];
end
df_ratings2 = table(repmat(noms2,size(exp2_data.rating_data,1),1),clr2,lrg2,ones(size(clr2,1),1));

% Reconcile items used in exp1 and exp2
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
df_keep = df_ratings(ismember(df_ratings.Var1,df_ratings2.Var1),:);
df_keep2 = df_ratings2(ismember(df_ratings2.Var1,df_ratings.Var1),:);
sum(~ismember(df_keep.Var1,df_keep2.Var1))
sum(~ismember(df_keep2.Var1,df_keep.Var1))
df_keep2.Properties.VariableNames = {'Var1','clr','lrg','Var4'};
df_final = [df_keep;df_keep2];
cd(dr.rpro)
writetable(df_final,'df_permutation_ratings.csv')

%% Correlations for response to reviewers analysis
% imagesc(corr(rating_color{:,2:end}'))
% idx = item_ind.cat1 == 1 | item_ind.cat1 == 3;

% Remove clr/lrg prefix
for i_row = 1:size(rating_large,1)
    rating_color.OriginalVariableNames{i_row} = rating_color.OriginalVariableNames{i_row}(5:end);
    rating_large.OriginalVariableNames{i_row} = rating_large.OriginalVariableNames{i_row}(5:end);
end

% Grab arranged names
color_noms_arranged = noms_arranged.noms.color;
large_noms_arranged = noms_arranged.noms.large;
all_noms_arranged = [color_noms_arranged; large_noms_arranged];

% Extract ratings across items
new_rating_struct = nan(size(all_noms_arranged,1),size(rating_color,2)-1);
for i_tem = 1:size(all_noms_arranged,1)
    tmp_nom = cell2mat(all_noms_arranged{i_tem}); % target name
    for i_rat = 1:size(rating_color,1)
        tmp_nom_check_color = rating_color.OriginalVariableNames{i_rat}; % check name
        tmp_nom_check_large = rating_large.OriginalVariableNames{i_rat}; % check name
        if strcmp(tmp_nom, tmp_nom_check_color) || strcmp(tmp_nom, tmp_nom_check_large)
            new_rating_struct(i_tem,:) = table2array(rating_color(i_rat,2:end));
        end
    end
end

% Add correlation and distance into df
correlation_matrix_ratings_color = corr(new_rating_struct(1:42,:)');
correlation_matrix_ratings_large = corr(new_rating_struct(43:end,:)');

% Check correlation matrix
figure
tiledlayout(2,2,"TileSpacing","compact")
nexttile(1)
imagesc(correlation_matrix_ratings_color)
colorbar
title('correlation - color')
nexttile(2)
imagesc(correlation_matrix_ratings_large)
colorbar
title('correlation - large')
nexttile(3)
imagesc(1-correlation_matrix_ratings_color)
colorbar
title('distance - color')
nexttile(4)
imagesc(1-correlation_matrix_ratings_large)
colorbar
title('distance - large')

% Vectorize matrices before computing distance
correlation_matrix_ratings_color(logical(eye(size(correlation_matrix_ratings_color)))) = 0;
correlation_vector_ratings_color = squareform(correlation_matrix_ratings_color);
distance_vector_ratings_color = 1-correlation_vector_ratings_color;
correlation_vector_ratings_color = atanh(correlation_vector_ratings_color);
distance_vector_ratings_color = atanh(distance_vector_ratings_color);
correlation_matrix_ratings_large(logical(eye(size(correlation_matrix_ratings_large)))) = 0;
correlation_vector_ratings_large = squareform(correlation_matrix_ratings_large);
distance_vector_ratings_large = 1-correlation_vector_ratings_large;
correlation_vector_ratings_large = atanh(correlation_vector_ratings_large);
distance_vector_ratings_large = atanh(distance_vector_ratings_large);

%% Save df
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
df = readtable('df_regression_exp2.csv');
df.Var11 = [correlation_vector_ratings_color';correlation_vector_ratings_large'];
df.Var12 = [distance_vector_ratings_color';distance_vector_ratings_large'];
cd(dr.rpro)
writetable(df,'df_regression_exp2.csv');