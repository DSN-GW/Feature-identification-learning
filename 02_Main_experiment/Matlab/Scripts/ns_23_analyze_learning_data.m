%% Analyze learning data
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
load('nsDataSort_after_exclusions.mat')

% Color info for visualizations
clrs = [     0, 0.4470, 0.7410; 
        0.8500, 0.3250, 0.0980;
        0.9290, 0.6940, 0.1250;
        0.4940, 0.1840, 0.5560;
        0.4660, 0.6740, 0.1880;
        0.3010, 0.7450, 0.9330;
        0.6350, 0.0780, 0.1840];
  
    
%% First, update sorting
for i_sub = 1:size(nsDataSort.data_all_scat,2)
    nsDataSort.data_all_scat{i_sub} = sortrows(nsDataSort.data_all_scat{i_sub},[6,4]);
end

%% Gather up PEs
PEs.color = [];
PEs.large = [];
gender = [];
for i_sub = 1:size(nsDataSort.data_all_scat,2)
    sub_data = nsDataSort.data_all_scat{i_sub};
    gender = [gender,sub_data.gender(1)];
    idx = sub_data.profile == 1;
    PEs.color = [PEs.color,sub_data.PEabsolute(idx)];
    PEs.large = [PEs.large,sub_data.PEabsolute(~idx)];
end


%% Compute means and standard errors
means.color = nanmean(PEs.color,2);
means.large = nanmean(PEs.large,2);
serrs.color = nanstd(PEs.color,[],2)/sqrt(size(PEs.color,2));
serrs.large = nanstd(PEs.large,[],2)/sqrt(size(PEs.large,2));

%% Visualize curves
ntrials = [1:size(PEs.color,1)]';
yup = [0,5];
xup = [0,43];
figpos = [10 10 1200 500];
figure('Position',figpos)
subplot(1,2,1)
boundedline(ntrials,means.color,serrs.color,'cmap',clrs(1,:))
[r,p] = corrcoef(means.color,ntrials);
ply = polyfit(ntrials,means.color,1);
ply = polyval(ply,ntrials);
hold on
h = plot(ntrials,ply,'LineWidth',2);
xlim(xup)
ylim(yup)
ylabel('absolute prediction error')
xlabel('trial number')
title('colorful')
legend(h,{[' r = ',num2str(r(1,2),2), ', p = ', num2str(p(1,2),2)]},'box','off')
set(gca,'FontSize',30)
subplot(1,2,2)
boundedline(ntrials,means.large,serrs.large,'cmap',clrs(2,:))
[r,p] = corrcoef(means.large,ntrials);
ply = polyfit(ntrials,means.large,1);
ply = polyval(ply,ntrials);
hold on
h = plot(ntrials,ply,'Color',clrs(2,:),'LineWidth',2);
xlim(xup)
ylim(yup)
title('large')
xlabel('trial number')
legend(h,{[' r = ', num2str(r(1,2),2), ', p = ', num2str(p(1,2),2)]},'box','off')
set(gca,'FontSize',30)

%% Gender comparison
idx1 = gender==1;
idx2 = gender==2;
idx3 = gender==3;

% Compute means
gen1_mean.color = nanmean(PEs.color(:,idx1),2);
gen2_mean.color = nanmean(PEs.color(:,idx2),2);
gen1_mean.large = nanmean(PEs.large(:,idx1),2);
gen2_mean.large = nanmean(PEs.large(:,idx2),2);
gen3_mean.color = nanmean(PEs.color(:,idx3),2);
gen3_mean.large = nanmean(PEs.large(:,idx3),2);

% Visualize
fig = figure;
subplot(1,2,1)
title('colorful')
hold on
plot(1:42,gen1_mean.color)
plot(1:42,gen2_mean.color)
plot(1:42,gen3_mean.color)
ylim([0,6])
xlim([0,43])
box on
legend({'M','W',"NB"})
subplot(1,2,2)
title('large')
hold on
plot(1:42,gen1_mean.large)
plot(1:42,gen2_mean.large)
plot(1:42,gen3_mean.large)
ylim([0,6])
box on
xlim([0,43])
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'prediction error');
xlabel(han,'trial number');

% Compute means
gen1_subject_mean.color = nanmean(PEs.color(:,idx1),1);
gen1_subject_mean.large = nanmean(PEs.large(:,idx1),1);
gen2_subject_mean.color = nanmean(PEs.color(:,idx2),1);
gen2_subject_mean.large = nanmean(PEs.color(:,idx2),1);
gen3_subject_mean.color = nanmean(PEs.color(:,idx3),1);
gen3_subject_mean.large = nanmean(PEs.large(:,idx3),1);

% Create table
tbl = table([gender';gender'],[repmat(0,length(gender),1);repmat(1,length(gender),1)], ...
             [gen1_subject_mean.color';gen2_subject_mean.color';gen3_subject_mean.color'; ...
              gen1_subject_mean.large';gen2_subject_mean.large';gen3_subject_mean.large']);
tbl.Properties.VariableNames = {'gender','condition','mean'};
tbl.gender = categorical(tbl.gender);
tbl.condition = categorical(tbl.condition);

% Comparison
mdl = fitlm(tbl,'mean ~ gender * condition')
anova(mdl)

% Save table
cd([dr.top,dr.exp,'/Data'])
writetable(tbl,'df_gender_PEs.csv')

% Trial level
tbl2 = table([repelem(1:82,42)';repelem(1:82,42)'], ...
             [repelem(gender,42)';repelem(gender,42)'], ...
             [repmat(0,3444,1);repmat(1,3444,1)], ...
             [PEs.color(:);PEs.large(:)], ...
             repmat(ntrials,164,1))
tbl2.Properties.VariableNames = {'subid','gender','condition','error','trial'};
tbl2.gender = categorical(tbl2.gender);
tbl2.condition = categorical(tbl2.condition);
mdl = fitlm(tbl2,'error ~ gender * condition')
anova(mdl)
cd([dr.top,dr.exp,'/Data'])
writetable(tbl2,'df_gender_PEs_trials.csv')

%% Checking how many subjects had significant PE reduction over time
clr_res = [];
lrg_res = [];
figure
for i_sub = 1:size(PEs.color,2)
    [r,p] = corrcoef(PEs.color(:,i_sub),ntrials);
    clr_res(i_sub,1) = r(1,2);
    clr_res(i_sub,2) = p(1,2);

    subplot(2,2,1)
    title('colorful')
    plot(ntrials,PEs.color(:,i_sub))
    hold on

    [r,p] = corrcoef(PEs.large(:,i_sub),ntrials);
    lrg_res(i_sub,1) = r(1,2);
    lrg_res(i_sub,2) = p(1,2);

    subplot(2,2,2)
    title('large')
    plot(ntrials,PEs.large(:,i_sub))
    hold on
end
sum(clr_res(:,2)<0.05 & clr_res(:,1)<0)
sum(lrg_res(:,2)<0.05 & lrg_res(:,1)<0)
sgtitle('subject-level curves')
subplot(2,2,3)
plot(ntrials,PEs.large(:,1))
ylim([0,10])
ply = polyfit(ntrials,PEs.large(:,1),1);
ply = polyval(ply,ntrials);
hold on
h = plot(ntrials,ply,'k--','LineWidth',2);
[r,p] = corrcoef(ntrials,PEs.large(:,1));
title('example subject')
legend('subject','fit = -0.12, p = 0.45')
subplot(2,2,4)
plot(ntrials,nanmean(PEs.color,2),'r');
title('colorful')
hold on
ply = polyfit(ntrials,nanmean(PEs.color,2),1);
ply = polyval(ply,ntrials);
hold on
h = plot(ntrials,ply,'r--','LineWidth',2);
[r,p] = corrcoef(ntrials,nanmean(PEs.color,2));
plot(ntrials,nanmean(PEs.large,2),'b')
ply = polyfit(ntrials,nanmean(PEs.large,2),1);
ply = polyval(ply,ntrials);
hold on
h = plot(ntrials,ply,'b--','LineWidth',2);
[r,p] = corrcoef(ntrials,nanmean(PEs.large,2));
ylim([0,10])
title('averaged curves')
legend({'colorful','fit1 = -0.61, p = 0.00','large','fit2 = -0.84, p = 0.00'})