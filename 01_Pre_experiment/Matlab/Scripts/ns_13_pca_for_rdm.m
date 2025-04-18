%% Run principal components analysis on RDM data
clear
close all

% Set directories
matpath = matlab.desktop.editor.getActiveFilename;
idx     = strfind(matpath,'\');
dr.top  = matpath(1:idx(end-2)-1);
dr.exp  = '/Matlab';
dr.data = '/Data';
dr.rpro = '/R';
rpath   = [dr.top,dr.rpro,'/Data'];
addpath(genpath(dr.top));
cd([dr.top,dr.exp,dr.data]);

% Read in data and stimulus info
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
load('allrdms_exp1.mat') % Data
load('rdm_sqmat.mat')    % Data
load('item_ind_all.mat') % Stimulus info

%% First steps
% Initialize a matrix with average dissimilarity for each item with all the
% others per participant
rdm_mean(1:120,1:size(rdm_mats,1)) = nan;

% Next, z-score dissimilarity values within participants
for i_sub = 1:size(rdm_mats,1)
    rdm = normalize(rdm_mats{i_sub,1}.data);
    rdm_mean(:,i_sub) = mean(rdm,2,'omitnan');
end


% Exclude participants with extensive missings
% Checking proportion of NaN from the raw RDM vectors for exclusions
nanper = [];
for i_sub = 1:size(allrdms,2)
    rdm = allrdms(:,i_sub);
    nanper(i_sub,1) = sum(isnan(rdm))/length(rdm);
end
histogram(nanper)
% mean(nanper)
% std(nanper)
% mean(nanper) + (std(nanper)*2)
% sum(nanper>0.5346)

% Average proportion of missing item-pair dissimilarities is 25% (SD = 16%)
% Using a cutoff of mean + 2SD results in a loss of 8 subjects
exclude_sub = nanper'>(mean(nanper)+(2*std(nanper)));

% Identify subjects who should be removed
% exclude_ids = data.meadows.meta(exclude_sub,10);
% save('exclude_ids.mat','exclude_ids')

% [row,sub]=find(isnan(rdm_mean));
% exclude_sub = unique(sub); % excluded subs
rdm_mean(:,exclude_sub)=[];

% Next, flip the matrix for PCA:columns vars and rows participants
rdm_PCA=rdm_mean';

%% PCA
[PCA.coeff, PCA.score, PCA.latent, PCA.tsquared, PCA.explained, PCA.mu] = pca(rdm_PCA);

% Calculate eigenvalues/vectors from covariance matrix
covariance_matrix = cov(rdm_PCA);
[V, D] = eig(covariance_matrix);

% How many components to retain?
cumexp = cumsum(PCA.explained);
ncomps_k  = sum(diag(D)>3); % Kaiser criterion
ncomps_90 = sum(cumexp<90)+1;

%% Reconstruct mat from components
recon1 = corrcoef(PCA.score(:,1:1)*PCA.coeff(:,1:1)'); % n_factor is 10 here
recon2 = corrcoef(PCA.score(:,1:2)*PCA.coeff(:,1:2)'); % n_factor is 10 here
recon3 = corrcoef(PCA.score(:,1:3)*PCA.coeff(:,1:3)');
recon12 = corrcoef(PCA.score(:,1:12)*PCA.coeff(:,1:12)');
recon_all = corrcoef(PCA.score*PCA.coeff');

% Fix underscores for axis labels
noms = item_ind.filename;
for i_tem = 1:size(noms,1)
    noms(i_tem) = strrep(noms{i_tem},'_','');
end
ind = zeros(119,1);
ind(1:5:119) = 1;
ind = [ind(2:end);0;1];

% Visualize
figpos = [20 20 1600 600];
figure('Position',figpos)
subplot(1,2,1)
imagesc(recon1);colorbar
xticks(0:5:120)
yticks(0:5:120)
xticklabels(noms(ind==1))
yticklabels(noms(ind==1))
xtickangle(90)
title({'Reconstruction using';'first principal component'})
set(gca,'FontSize',15)
subplot(1,2,2)
imagesc(recon2);colorbar
xticks(0:5:120)
yticks(0:5:120)
xticklabels(noms(ind==1))
yticklabels(noms(ind==1))
xtickangle(90)
title({'Reconstruction using';'first two principal components'})
set(gca,'FontSize',15)

%% Build a category mat that can be tested against components
cat=(item_ind.cat1);

% Primary category
for i=1:size(unique(cat))
    sumcat(i)= sum(cat==i);
end
for i=1:size(unique(cat))
    sumcat(2,i)=sumcat(1,i)*(sumcat(1,i)-1)/2;
end
catvec=zeros(sum(120*(120-1)/2),1);
cat_mat=squareform(catvec);
catvec1=ones(sumcat(1,1),sumcat(1,1));
catvec2=ones(sumcat(1,2),sumcat(1,2));
catvec3=ones(sumcat(1,3),sumcat(1,3));
cat_mat(1:sumcat(1,1),1:sumcat(1,1))=catvec1;
cat_mat(sumcat(1,1)+1:sum(sumcat(1,1:2)),sumcat(1,1)+1:sum(sumcat(1,1:2)))=catvec2;
cat_mat(sum(sumcat(1,1:2))+1:end,sum(sumcat(1,1:2))+1:end)=catvec3;

% Subcategory
cat_mat_2 = zeros(120,120);
cat2 = item_ind.cat2;
for i=1:size(unique(cat2))
    sumcat2(i) = sum(cat2==i);
end
cat_mat_2(1:sumcat2(i),1:sumcat2(i)) = ones(sumcat2(i),sumcat2(i));
for i=1:size(sumcat2,2)
    cat_grid = ones(sumcat2(i),sumcat2(i));     %create subcat grid
    stt = 1 + sum(sumcat2(1:i)) - sumcat2(i);   %start index
    stp = stt + sumcat2(i) - 1;                 %stop index
    cat_mat_2(stt:stp,stt:stp) = cat_grid;      %place grid in cat_mat_2
end
imagesc(cat_mat_2); colorbar

% Turn cat_mat to one vector
cat_low=tril(cat_mat);
cat_array=cat_low(:);
cat2_low = tril(cat_mat_2);
cat2_array = cat2_low(:);

% Vectorize
recon2_mat=tril(recon2);
recon2_array=recon2_mat(:);

%% Setting up df for regression modeling of reconstructions
reconstructions = [];
for i_tem = 1:size(PCA.score,2)
    rc = corrcoef(PCA.score(:,1:i_tem)*PCA.coeff(:,1:i_tem)');
    %     imagesc(rc)
    %     pause()
    %     close
    rc = tril(rc);
    rc = rc(:);
    reconstructions(:,i_tem) = rc;
end
reconstructions = [cat_array,cat2_array,reconstructions];

% Simple regression df
m_rdm = nanmean(allrdms,2);
m_rdm = squareform(m_rdm);
m_rdm_norm = normalize(m_rdm);

% Visualize
figure
imagesc(m_rdm_norm); colorbar

% Normalize RDMs
n_rdm = [];
for i_sub = 1:size(allrdms,2)
    n_rdm(:,i_sub) = normalize(allrdms(:,i_sub));
end
n_rdm = squareform(nanmean(n_rdm,2));

% Visualize
figure
imagesc(n_rdm); colorbar
m_rdm_low = tril(squareform(m_rdm));
m_rdm_array = m_rdm_low(:);
n_rdm_low = tril(n_rdm);
n_rdm_array = n_rdm_low(:);
figure
imagesc(m_rdm - n_rdm); colorbar
df_regression_simple = [n_rdm_array,cat_array,cat2_array];

%% Save mats and PCA mats
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
writematrix(reconstructions,[rpath,'/df_regression_exp1.csv'])
writematrix(df_regression_simple,[rpath,'/df_regression_simple_exp1.csv'])
save(fullfile([dr.top,dr.exp,'/data/', 'PCA_mats']),'PCA', 'rdm_mean', ...
    'rdm_PCA','exclude_sub');

%% More visualizations
figpos = [10 10 700 600];
figure('Position',figpos)
bar([1:size(D,1)],sort(diag(D),'descend'));
set(gca,'FontSize',15)
title('eigenvalues')

% Imagesc for recon
figpos = [10 10 1100 900];
figure('Position',figpos)
subplot(2,2,1)
imagesc(recon1); colorbar
ylabel('item number')
title('reconstruction on 1 PC')
set(gca,'FontSize',15)
subplot(2,2,2)
imagesc(recon3); colorbar
title('reconstruction on 3 PCs')
set(gca,'FontSize',15)
subplot(2,2,3)
imagesc(recon12); colorbar
xlabel('item number')
ylabel('item number')
title('reconstruction on 12 PCs')
set(gca,'FontSize',15)
subplot(2,2,4)
% imagesc(recon_all); colorbar
% xlabel('item number')
% title('reconstruction on all PCs')
% set(gca,'FontSize',15)
bar([1:120],cumexp)
xlabel('principal component')
ylim([0,100])
ylabel('variance explained (%)')
title('cumulative distribution of variance explained')
set(gca,'FontSize',15)

%% More visualizations
% Feature loadings
figpos = [10 10 900 1300];
figure('Position',figpos)
h = imagesc(PCA.coeff); colorbar
xlabel('principal componenet')
yticks(1:120)
yticklabels(item_ind.filename)
title('feature loadings')

% First 5 components
figpos = [10 10 900 1300];
figure('Position',figpos)
h = imagesc(PCA.coeff(:,1:5)); colorbar
xlabel('principal componenet')
xticks(1:1:5)
yticks(1:120)
yticklabels(item_ind.filename)
title(['feature loadings (first five components, ', ...
    num2str(cumexp(5),2),'% variance explained)'])

% Scree plot and cumulative variance
figpos = [10 10 1200 500];
figure('Position',figpos)
subplot(1,2,1)
bar([1:120],PCA.explained)
xlabel('principal component')
ylabel('variance explained (%)')
ylim([-1,100])
title('scree plot')
set(gca,'FontSize',15)
subplot(1,2,2)
bar([1:120],cumexp)
xlabel('principal component')
ylim([-1,100])
title('cumulative')
set(gca,'FontSize',15)

% Set colors
subclrs = [[204   0   0]/255;
    [204   0   0]/255;
    [204   0   0]/255;
    [204   0   0]/255;
    [204 204   0]/255;
    [204 204   0]/255;
    [204 204   0]/255;
    [204 204   0]/255;
    [  0 102 204]/255;
    [  0 102 204]/255;
    [  0 102 204]/255;
    [  0 102 204]/255];

% Data projection onto first 3 components
figpos = [10 10 700 700];
figure('Position',figpos)
for i_tem = 1:size(item_ind,1)
    if item_ind.cat1(i_tem) == 1
        h1 = scatter3(PCA.coeff(i_tem,1),PCA.coeff(i_tem,2),PCA.coeff(i_tem,3), ...
            'MarkerFaceColor',subclrs(1,:),'MarkerEdgeColor',subclrs(1,:));
        hold on
    elseif item_ind.cat1(i_tem) == 2
        h2 = scatter3(PCA.coeff(i_tem,1),PCA.coeff(i_tem,2),PCA.coeff(i_tem,3), ...
            'MarkerFaceColor',subclrs(5,:),'MarkerEdgeColor',subclrs(5,:));
        hold on
    elseif item_ind.cat1(i_tem) == 3
        h3 = scatter3(PCA.coeff(i_tem,1),PCA.coeff(i_tem,2),PCA.coeff(i_tem,3), ...
            'MarkerFaceColor',subclrs(9,:),'MarkerEdgeColor',subclrs(9,:));
        hold on
    end
end
xlabel('1st component')
ylabel('2nd component')
zlabel('3rd component')
title('data projected onto first three components')
legend([h1,h2,h3],{'activities','fashion','foods'},'Location','northwest')
set(gca,'FontSize',15)

% Just two components
figpos = [10 10 500 400];
figure('Position',figpos)
for i_tem = 1:size(item_ind,1)
    if item_ind.cat1(i_tem) == 1
        h1 = scatter(PCA.coeff(i_tem,1),PCA.coeff(i_tem,2), ...
            'MarkerFaceColor',subclrs(1,:),'MarkerEdgeColor',subclrs(1,:));
        hold on
    elseif item_ind.cat1(i_tem) == 2
        h2 = scatter(PCA.coeff(i_tem,1),PCA.coeff(i_tem,2), ...
            'MarkerFaceColor',subclrs(5,:),'MarkerEdgeColor',subclrs(5,:));
        hold on
    elseif item_ind.cat1(i_tem) == 3
        h3 = scatter(PCA.coeff(i_tem,1),PCA.coeff(i_tem,2), ...
            'MarkerFaceColor',subclrs(9,:),'MarkerEdgeColor',subclrs(9,:));
        hold on
    end
end
xlabel('1st component')
ylabel('2nd component')
title('data projected onto first two components')
legend([h1,h2,h3],{'activities','fashion','foods'},'Location','northeast')
set(gca,'FontSize',15)

%% In response to reviewers
% Only frist two components, with only color/large items, for modeling
% Note that pre-experiment: exp1 and main experiment: exp2 in naming 

% Get arranged stimulus names
noms_arranged = load('arranged_names_exp2.mat');
color_noms_arranged = noms_arranged.noms.color;
large_noms_arranged = noms_arranged.noms.large;
all_noms_arranged = [color_noms_arranged; large_noms_arranged];
idx_color = zeros(120,1);
idx_large = zeros(120,1);

% Find index of items used in exp2
% Color
for i_tem = 1:size(item_ind,1)
    itm1 = item_ind.filename{i_tem};
    for i_row = 1:size(color_noms_arranged,1)
        itm2 = cell2mat(color_noms_arranged{i_row});
        if strcmp(itm1,itm2)
            idx_color(i_tem) = 1;
        end
    end
end
% Large
for i_tem = 1:size(item_ind,1)
    itm1 = item_ind.filename{i_tem};
    for i_row = 1:size(large_noms_arranged,1)
        itm2 = cell2mat(large_noms_arranged{i_row});
        if strcmp(itm1,itm2)
            idx_large(i_tem) = 1;
        end
    end
end

% From recon1 and 2, keep only items used in exp2
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
figure
tiledlayout(1,2,'TileSpacing','compact')
nexttile(1)
imagesc(recon1)
nexttile(2)
imagesc(recon2) % originals

% Logicals
idx_color = logical(idx_color);
idx_large = logical(idx_large);

% Extracted reconstructions
figure
tiledlayout(2,2,'TileSpacing','compact')
nexttile(1)
imagesc(recon1(idx_color,idx_color))
title('color')
ylabel('recon1')
nexttile(2)
imagesc(recon1(idx_large,idx_large))
title('large')
nexttile(3)
ylabel('recon2')
imagesc(recon2(idx_color,idx_color))
nexttile(4)
imagesc(recon2(idx_large,idx_large))

% unfold and export for modeling
% example:
% correlation_matrix_ratings_color(logical(eye(size(correlation_matrix_ratings_color)))) = 0;
clr_recon1 = recon1(idx_color,idx_color);
lrg_recon1 = recon1(idx_large,idx_large);
clr_recon2 = recon2(idx_color,idx_color);
lrg_recon2 = recon2(idx_large,idx_large);
clr_recon1(logical(eye(size(clr_recon1)))) = 0;
lrg_recon1(logical(eye(size(lrg_recon1)))) = 0;
clr_recon2(logical(eye(size(clr_recon2)))) = 0;
lrg_recon2(logical(eye(size(lrg_recon2)))) = 0;
extracted_recon1 = [squareform(clr_recon1)';squareform(lrg_recon1)'];
extracted_recon2 = [squareform(clr_recon2)';squareform(lrg_recon2)'];
cd([dr.top,dr.exp,dr.data]);
writematrix([extracted_recon1,extracted_recon2],'extracted_reconstructions.csv')