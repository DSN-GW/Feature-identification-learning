
%% family comparison script
%% Feature identification learning 
%% Computational modeling analyses for Real Data sets 
% #########################################################################
cat = load("./02_Main_experiment_computational_modeling/Data/data_model_fit_all13_v1cat.mat");
subcat = load("./02_Main_experiment_computational_modeling/Data/data_model_fit_all13_v1sub_cat.mat");

figure
bar([cat.bicMin',subcat.bicMin'])

%% add SPM in the path
addpath(genpath('./spm12'));
%% family comparison

re_arra = [3,4,5]; % compare main models with cat versus subcat resets
model = [];
lme = [ -cat.bicMatrix(:,re_arra), -subcat.bicMatrix(:,re_arra) ];
family.infer = 'RFX';
family.partition = [ ones(1,3), ones(1,3)+1 ];
family.names = {'fam1','fam2'};
[family, model] = spm_compare_families( lme, family );

figure 
bar(family.xp)

%% self vs mean RPs in main models
re_arra1 = [1,4,7];
re_arra2 = [2,5,8];
model = [];
lme = [ -cat.bicMatrix(:,re_arra1), -cat.bicMatrix(:,re_arra2) ];
family.infer = 'RFX';
family.partition = [ ones(1,3), ones(1,3)+1 ];
family.names = {'self','mean'};
[family, model] = spm_compare_families( lme, family );

figure
bar(family.xp)
