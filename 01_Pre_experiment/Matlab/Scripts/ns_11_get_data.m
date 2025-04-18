%% Read in raw data files and saves as data structure
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

% Read in list of IDs to exclude based on > M+2SD missing data
load('exclude_ids.mat')

%%  Read in raw data files from Prolific
prolific1 = readtable('prolific_export_606cb9026b4fe2da6318400e.csv'); % sample1
prolific2 = readtable('prolific_export_60c26f29978f91950fde551a.csv'); % sample2

% Find missing
empty_prolific1 = cellfun(@isempty,prolific1.participant_id);
empty_prolific2 = cellfun(@isempty,prolific2.participant_id);

% Check that none are missing
% sum([empty_prolific1;empty_prolific2]) == 0

% Set data structure
data.prolific.data = [prolific1(:,[1:10,12:19]);prolific2(:,[1:10,12:19])];
clear prolific1 prolific2

%% Read in raw data files from Redcap
redcap1 = readtable('AdultNonsocial_DATA_2021-04-08_1238.csv');  %sample1
redcap2 = readtable('AdultNonsocial2_DATA_2022-01-07_1507.csv'); %sample2

% Check if rows empty
empty_redcap1 = cellfun(@isempty,redcap1.subid);
empty_redcap2 = cellfun(@isempty,redcap2.subid);
redcap1 = redcap1(empty_redcap1==0,:);
redcap2 = redcap2(empty_redcap2==0,:);

%% Select common variables across the 2 redcap studies
idx1 = ismember(redcap1.Properties.VariableNames,redcap2.Properties.VariableNames); % reconciling
idx2 = ismember(redcap2.Properties.VariableNames,redcap1.Properties.VariableNames); % variable varnames
data.redcap.data = [redcap1(:,idx1);redcap2(:,idx2)];
clear redcap1 redcap2 idx1 idx2

%% Read in raw data files from Meadows
data.meadows.data = [struct2cell(jsondecode(fileread('Meadows_redcap_testing_02_v_v4_tree.json'))); ...
    struct2cell(jsondecode(fileread('Meadows_redcap_testing_02_v_v5_tree.json'))); ...
    struct2cell(jsondecode(fileread('Meadows_redcap_testing_03_v_v2_tree.json')))];
data.meadows.meta = [readtable('Meadows_redcap_testing_02_v_v4_metadata.csv'); ...
    readtable('Meadows_redcap_testing_02_v_v5_metadata.csv'); ...
    readtable('Meadows_redcap_testing_03_v_v2_metadata.csv')];

% Check no missings
empty_meadows = cellfun(@isempty,data.meadows.meta.subid);
% find(empty_meadows~=0)  % row 79 needs to be deleted
% sum(empty_meadows) == 0 % check no missings

%% Delete empty vars
clear empty*;

%% Keep only completed subids
data.prolific.data = data.prolific.data(strncmp(data.prolific.data{:,3},'A',1),:);
data.redcap.data   = data.redcap.data(data.redcap.data{:,391}==2,:);
data.meadows.meta  = data.meadows.meta(strcmp(data.meadows.meta{:,4},'finished'),:);

%% Adding missing prolific ids for meadows version 5
% get only v5 file
x = struct2cell(jsondecode(fileread('Meadows_redcap_testing_02_v_v5_tree.json')));

% prolific ids are stored in structure under task 1
for i_dee = 1:length(x)
    profid{i_dee,1} = x{i_dee}.tasks{1}.prolific_id;
end
clear i_dee

% add missing ids into full structure
data.meadows.meta{data.meadows.meta.version==5,10} = profid;
clear x profid

%%  check for duplicate subids
% Prolific
[D,~,X] = unique(data.prolific.data.participant_id(:));
Y = hist(X,unique(X));
sum(Y) == length(D); % TRUE means no duplicates
clear X Y D

% Redcap
[D,~,X] = unique(data.redcap.data.subid(:));
Y = hist(X,unique(X));
sum(Y) == length(D);
dupe1 = D(Y~=1); % One duplicate
data.redcap.data(strcmp(data.redcap.data.subid,dupe1),:);
for i_sub=1:size(data.redcap.data,1)
    dup_rows(i_sub)=sum(strcmp(data.redcap.data.subid{i_sub},dupe1));
end
  
% Delete duplicate rows
data.redcap.data(find(dup_rows~=0,1,'last'),:)=[];
clear X Y D dupe1 

% Meadows
[D,~,X] = unique(data.meadows.meta.subid(:));
Y = hist(X,unique(X));
sum(Y) == length(D);
clear X Y D

%% remove ids with too much missing data (M+2SD)
%  these subjects were identified in the PCA script
%  removing now and saving structure with exclusions
idx = ismember(data.redcap.data(:,4),exclude_ids);
data.redcap.data = data.redcap.data(~idx,:);
idx = ismember(data.meadows.meta(:,10),exclude_ids);
data.meadows.meta = data.meadows.meta(~idx,:);
data.meadows.data = data.meadows.data(~idx);
exclude_ids.Properties.VariableNames = {'participant_id'}; %need to match variable name with prolific
idx = ismember(data.prolific.data(:,2),exclude_ids);
data.prolific.data = data.prolific.data(~idx,:);

%% save data
% Note that pre-experiment: exp1 and main experiment: exp2 in naming
save(fullfile([dr.top,dr.exp,'/Data/', 'full_data_struct_exp1.mat']),'data');