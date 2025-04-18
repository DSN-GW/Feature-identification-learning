%% Read in raw data files and saves as data structure
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

% Stimulus info
load('item_ind_all.mat') 

%% read in raw data - sample 1
% Prolific
s1.prolific.data = readtable('sample1_prolific_export_6140d2bf320df68dc96327b1.csv');

% Redcap
s1.redcap.data = readtable('sample1_AdultNonsocialLearni_DATA_2021-09-15_1140.csv');

% Psiturk
s1.psiturk.data = readtable('sample1_trialdata.csv');

% Meadows
% reads in four files, metadata and actual data (ids) for color and large
% also specifies which batch of participants to select from. meadows generates
% a unique batch for each round of data collection. for this first sample, 
% we want only batch 'h'. 
s1.meadows.color.meta = readtable('sample1_Meadows_nonsocial_1_v_v1_metadata.csv');
s1.meadows.color.ids = jsondecode(fileread('sample1_Meadows_nonsocial_1_v_v1_tree.json'));
s1.meadows.large.meta = readtable('sample1_Meadows_nonsocial_2_v_v1_metadata.csv');
s1.meadows.large.ids = jsondecode(fileread('sample1_Meadows_nonsocial_2_v_v1_tree.json'));
s1.batch = 'h'; %% what is this?

%% Remove unneeded entries and sort for prolific, redcap, and psiturk
% in prolific, participants who have completed the experiment are listed as
% either 'AWAITING REVIEW' or 'APPROVED', depending on when the file is
% downloaded. Here, we want to select only those 'AWATING REVIEW'. We then
% sort the data according to subid (participant_id) so that all files align
% to subid
s1.prolific.data = s1.prolific.data(strcmp(s1.prolific.data.status,'AWAITING REVIEW'),:);
s1.prolific.data = sortrows(s1.prolific.data,2);
s1.prolific.data.time_taken = s1.prolific.data.time_taken./60; %time conversion

% In redcap, participants who have completed the experiment have a value 2
% (opposed to 0) in column 396, participant_ratings_complete. We want to select
% only those participants who completed the final rating task. then sort.
s1.redcap.data = s1.redcap.data(s1.redcap.data{:,396}==2,:);
s1.redcap.data = sortrows(s1.redcap.data,4);

% In psiturk, the data are sorted prior to import, so no need for sorting
% here. However, we want to remove missing rows/nonresponses, then we want a
% list of all of the unique subids
s1.psiturk.data = rmmissing(s1.psiturk.data);
s1.psiturk.subids = unique(s1.psiturk.data.subid);

%% Remove unneeded entries and sort for meadows - color
% Here, we start with the color data. First, we select the correct batch
% identified on line 33, then we select from those participants only those
% who completed the task, i.e., 'finished' in column 4 (status)
s1.meadows.color.meta = s1.meadows.color.meta(strcmp(s1.meadows.color.meta{:,3},s1.batch),:);
s1.meadows.color.meta = s1.meadows.color.meta(strcmp(s1.meadows.color.meta{:,4},'finished'),:);

% Here, we need to reconcile subids which are stored differently in meadows
% than in the other structures. In short, meadows gives its own unique id to
% each participant in the form of 'adjective-animal'. These ids are linked
% to the master subids in the metadata structure. So, what is done here:
% first, we use ismember() to crosscheck meadows ids in the metadata and
% real data structures to find which ids we can keep and which to remove
ids_to_keep = s1.meadows.color.meta.name;
ids_to_keep = strrep(ids_to_keep,'-','_'); %need to replace hyphen for strcmp()
idx.color = ismember(fieldnames(s1.meadows.color.ids),string(ids_to_keep));
fns.color = fieldnames(s1.meadows.color.ids);
fns.color = fns.color(~idx.color); %names to remove
s1.meadows.color.ids = rmfield(s1.meadows.color.ids,fns.color);

% Next, we sort the id list and the metadata file to match meadows ids with 
% subids. Then the data structure (ids) is sorted and the correct subid is
% added to the second column, leaving us we the correct subid associated
% with the correct participant entry. Finally, sort data by subid
ids_to_keep = sortrows(ids_to_keep,1);
s1.meadows.color.meta = sortrows(s1.meadows.color.meta,1);
ids_to_keep = [ids_to_keep,s1.meadows.color.meta.workerId];
s1.meadows.color.ids = orderfields(s1.meadows.color.ids);
s1.meadows.color.ids = struct2cell(s1.meadows.color.ids);
s1.meadows.color.ids = [s1.meadows.color.ids,ids_to_keep(:,2)];
s1.meadows.color.ids = sortrows(s1.meadows.color.ids,2);

%% Remove unneeded entries and sort for meadows - large
% Now, we repeat this process for the large data
s1.meadows.large.meta = s1.meadows.large.meta(strcmp(s1.meadows.large.meta.batch,s1.batch),:);
s1.meadows.large.meta = s1.meadows.large.meta(strcmp(s1.meadows.large.meta.status,'finished'),:);
ids_to_keep = s1.meadows.large.meta.name;
ids_to_keep = strrep(ids_to_keep,'-','_');
idx.large = ismember(string(fieldnames(s1.meadows.large.ids)),string(ids_to_keep));
fns.large = fieldnames(s1.meadows.large.ids);
fns.large = fns.large(~idx.large);
s1.meadows.large.ids = rmfield(s1.meadows.large.ids,fns.large);
ids_to_keep = sortrows(ids_to_keep);
s1.meadows.large.meta = sortrows(s1.meadows.large.meta,1);
ids_to_keep = [ids_to_keep,s1.meadows.large.meta.workerId];
s1.meadows.large.ids = orderfields(s1.meadows.large.ids);
s1.meadows.large.ids = struct2cell(s1.meadows.large.ids);
s1.meadows.large.ids = [s1.meadows.large.ids,ids_to_keep(:,2)];
s1.meadows.large.ids = sortrows(s1.meadows.large.ids,2);

%% Reconcile ids
% Now, we need to crosscheck subids across each data source so that we know
% who completed all tasks.

% First, we check for color participants who completed large
idx.keep = ismember(s1.meadows.color.ids(:,2),s1.meadows.large.ids(:,2));
s1.meadows.color.meta = s1.meadows.color.meta(idx.keep,:);
s1.meadows.color.ids = s1.meadows.color.ids(idx.keep,:);
% Then large participants who completed color
idx.keep = ismember(s1.meadows.large.ids(:,2),s1.meadows.color.ids(:,2));
s1.meadows.large.meta = s1.meadows.large.meta(idx.keep,:);
s1.meadows.large.ids = s1.meadows.large.ids(idx.keep,:);

% Next, which meadows participants completed redcap
idx.keep = ismember(s1.meadows.color.ids(:,2),s1.redcap.data.subid);
s1.meadows.color.meta = s1.meadows.color.meta(idx.keep,:);
s1.meadows.color.ids = s1.meadows.color.ids(idx.keep,:);
s1.meadows.large.meta = s1.meadows.large.meta(idx.keep,:);
s1.meadows.large.ids = s1.meadows.large.ids(idx.keep,:);
% And which redcap participants completed meadows
idx.keep = ismember(s1.redcap.data.subid,s1.meadows.color.ids(:,2));
s1.redcap.data = s1.redcap.data(idx.keep,:);

% Next, which psiturk participants completed redcap
idx.keep = ismember(s1.psiturk.subids,s1.redcap.data.subid);
s1.psiturk.subids = s1.psiturk.subids(idx.keep);
s1.psiturk.data = s1.psiturk.data(ismember(s1.psiturk.data.subid,s1.psiturk.subids),:);
% And which redcap participants completed psiturk
idx.keep = ismember(s1.redcap.data.subid,s1.psiturk.subids);
s1.redcap.data = s1.redcap.data(idx.keep,:);
% And which meadows participants completed psiturk
idx.keep = ismember(s1.meadows.color.ids(:,2),s1.psiturk.subids);
s1.meadows.color.meta = s1.meadows.color.meta(idx.keep,:);
s1.meadows.color.ids = s1.meadows.color.ids(idx.keep,:);
s1.meadows.large.meta = s1.meadows.large.meta(idx.keep,:);
s1.meadows.large.ids = s1.meadows.large.ids(idx.keep,:);

% Finally, which prolific participants completed psiturk
idx.keep = ismember(s1.prolific.data.participant_id,s1.psiturk.subids);
s1.prolific.data = s1.prolific.data(idx.keep,:);
% And which psiturk participants completed prolific
idx.keep = ismember(s1.psiturk.subids,s1.prolific.data.participant_id);
s1.psiturk.subids = s1.psiturk.subids(idx.keep);
s1.psiturk.data = s1.psiturk.data(ismember(s1.psiturk.data.subid,s1.psiturk.subids),:);
% And which redcap participants completed prolific
idx.keep = ismember(s1.redcap.data.subid,s1.prolific.data.participant_id);
s1.redcap.data = s1.redcap.data(idx.keep,:);
% And which meadows participants completed prolific
idx.keep = ismember(s1.meadows.color.ids(:,2),s1.prolific.data.participant_id);
s1.meadows.color.meta = s1.meadows.color.meta(idx.keep,:);
s1.meadows.color.ids = s1.meadows.color.ids(idx.keep,:);
s1.meadows.large.meta = s1.meadows.large.meta(idx.keep,:);
s1.meadows.large.ids = s1.meadows.large.ids(idx.keep,:);

% Now that we have the full, organized structure for our first sample, we
% repeat the process for our second sample. once done, we combine & clean

%% Read in raw data - sample 2
s2.prolific.data = readtable('sample2_prolific_export_6165a5744ac311682252f608.csv');
s2.redcap.data = readtable('sample2_AdultNonsocialLearni_DATA_2021-10-14_1337.csv');
s2.psiturk.data = readtable('sample2_trialdata.csv');
s2.meadows.color.meta = readtable('sample2_Meadows_nonsocial_1_v_v4_metadata.csv');
s2.meadows.color.ids = jsondecode(fileread('sample2_Meadows_nonsocial_1_v_v4_tree.json'));
s2.meadows.large.meta = readtable('sample2_Meadows_nonsocial_2_v_v4_metadata.csv');
s2.meadows.large.ids = jsondecode(fileread('sample2_Meadows_nonsocial_2_v_v4_tree.json'));
s2.batch = 'a'; %for sample 2 we want batch 'a' instead of 'h'

% Remove unneeded entries and sort
% note: psiturk sorted before import
s2.prolific.data = s2.prolific.data(strcmp(s2.prolific.data.status,'AWAITING REVIEW'),:); %finished?
s2.prolific.data = sortrows(s2.prolific.data,2);
s2.prolific.data.time_taken = s2.prolific.data.time_taken./60; %time conversion
s2.redcap.data = s2.redcap.data(s2.redcap.data{:,401}==2,:); %finished?
s2.redcap.data = sortrows(s2.redcap.data,4);
s2.psiturk.data = rmmissing(s2.psiturk.data);
s2.psiturk.subids = unique(s2.psiturk.data.subid);

% Meadows - color
s2.meadows.color.meta = s2.meadows.color.meta(strcmp(s2.meadows.color.meta{:,3},s2.batch),:); %which batch
s2.meadows.color.meta = s2.meadows.color.meta(strcmp(s2.meadows.color.meta{:,4},'finished'),:); %finished?
ids_to_keep = s2.meadows.color.meta.name;
ids_to_keep = strrep(ids_to_keep,'-','_');
idx.color = ismember(fieldnames(s2.meadows.color.ids),string(ids_to_keep));
fns.color = fieldnames(s2.meadows.color.ids);
fns.color = fns.color(~idx.color); %names to remove
s2.meadows.color.ids = rmfield(s2.meadows.color.ids,fns.color);
ids_to_keep = sortrows(ids_to_keep,1);
s2.meadows.color.meta = sortrows(s2.meadows.color.meta,1);
ids_to_keep = [ids_to_keep,s2.meadows.color.meta.workerId];
s2.meadows.color.ids = orderfields(s2.meadows.color.ids);
s2.meadows.color.ids = struct2cell(s2.meadows.color.ids);
s2.meadows.color.ids = [s2.meadows.color.ids,ids_to_keep(:,2)];
s2.meadows.color.ids = sortrows(s2.meadows.color.ids,2);
s2.meadows.color.ids = s2.meadows.color.ids(1:end-2,:);
s2.meadows.color.meta = sortrows(s2.meadows.color.meta,15);
s2.meadows.color.meta = s2.meadows.color.meta(1:end-2,:);

% Meadows - large
s2.meadows.large.meta = s2.meadows.large.meta(strcmp(s2.meadows.large.meta.batch,s2.batch),:);
s2.meadows.large.meta = s2.meadows.large.meta(strcmp(s2.meadows.large.meta.status,'finished'),:);
ids_to_keep = s2.meadows.large.meta.name;
ids_to_keep = strrep(ids_to_keep,'-','_');
idx.large = ismember(string(fieldnames(s2.meadows.large.ids)),string(ids_to_keep));
fns.large = fieldnames(s2.meadows.large.ids);
fns.large = fns.large(~idx.large);
s2.meadows.large.ids = rmfield(s2.meadows.large.ids,fns.large);
ids_to_keep = sortrows(ids_to_keep);
s2.meadows.large.meta = sortrows(s2.meadows.large.meta,1);
ids_to_keep = [ids_to_keep,s2.meadows.large.meta.workerId];
s2.meadows.large.ids = orderfields(s2.meadows.large.ids);
s2.meadows.large.ids = struct2cell(s2.meadows.large.ids);
s2.meadows.large.ids = [s2.meadows.large.ids,ids_to_keep(:,2)];
s2.meadows.large.ids = sortrows(s2.meadows.large.ids,2);
s2.meadows.large.ids = s2.meadows.large.ids(1:end-3,:);
s2.meadows.large.meta = sortrows(s2.meadows.large.meta,14);
s2.meadows.large.meta = s2.meadows.large.meta(1:end-3,:);

% Reconcile ids
% Within meadows
idx.keep = ismember(s2.meadows.color.ids(:,2),s2.meadows.large.ids(:,2));
s2.meadows.color.meta = s2.meadows.color.meta(idx.keep,:);
s2.meadows.color.ids = s2.meadows.color.ids(idx.keep,:);
idx.keep = ismember(s2.meadows.large.ids(:,2),s2.meadows.color.ids(:,2));
s2.meadows.large.meta = s2.meadows.large.meta(idx.keep,:);
s2.meadows.large.ids = s2.meadows.large.ids(idx.keep,:);

% Meadows and redcap
idx.keep = ismember(s2.meadows.color.ids(:,2),s2.redcap.data.subid);
s2.meadows.color.meta = s2.meadows.color.meta(idx.keep,:);
s2.meadows.color.ids = s2.meadows.color.ids(idx.keep,:);
s2.meadows.large.meta = s2.meadows.large.meta(idx.keep,:);
s2.meadows.large.ids = s2.meadows.large.ids(idx.keep,:);
idx.keep = ismember(s2.redcap.data.subid,s2.meadows.color.ids(:,2));
s2.redcap.data = s2.redcap.data(idx.keep,:);

% Now with psiturk
idx.keep = ismember(s2.psiturk.subids,s2.redcap.data.subid);
s2.psiturk.subids = s2.psiturk.subids(idx.keep);
s2.psiturk.data = s2.psiturk.data(ismember(s2.psiturk.data.subid,s2.psiturk.subids),:);
idx.keep = ismember(s2.redcap.data.subid,s2.psiturk.subids);
s2.redcap.data = s2.redcap.data(idx.keep,:);
idx.keep = ismember(s2.meadows.color.ids(:,2),s2.psiturk.subids);
s2.meadows.color.meta = s2.meadows.color.meta(idx.keep,:);
s2.meadows.color.ids = s2.meadows.color.ids(idx.keep,:);
s2.meadows.large.meta = s2.meadows.large.meta(idx.keep,:);
s2.meadows.large.ids = s2.meadows.large.ids(idx.keep,:);

% Now prolific
idx.keep = ismember(s2.prolific.data.participant_id,s2.psiturk.subids);
s2.prolific.data = s2.prolific.data(idx.keep,:);
idx.keep = ismember(s2.psiturk.subids,s2.prolific.data.participant_id);
s2.psiturk.subids = s2.psiturk.subids(idx.keep);
s2.psiturk.data = s2.psiturk.data(ismember(s2.psiturk.data.subid,s2.psiturk.subids),:);
idx.keep = ismember(s2.redcap.data.subid,s2.prolific.data.participant_id);
s2.redcap.data = s2.redcap.data(idx.keep,:);
idx.keep = ismember(s2.meadows.color.ids(:,2),s2.prolific.data.participant_id);
s2.meadows.color.meta = s2.meadows.color.meta(idx.keep,:);
s2.meadows.color.ids = s2.meadows.color.ids(idx.keep,:);
s2.meadows.large.meta = s2.meadows.large.meta(idx.keep,:);
s2.meadows.large.ids = s2.meadows.large.ids(idx.keep,:);

%% Combine samples 1 and 2 into overall data structure
data.prolific.data = [s1.prolific.data;s2.prolific.data];
data.redcap.data = [s1.redcap.data;s2.redcap.data];
% Psiturk note: due to modifcations to psiturk between sample collections, 
% columns from samples 1 and 2 do not perf. match, so we keep only matching 
% columns and we update column 4 in sample 1 from 'image_name' to 'filename'
s2.psiturk.data = [s2.psiturk.data(:,1),s2.psiturk.data(:,4:end)];
s1.psiturk.data.Properties.VariableNames{4} = 'filename';
data.psiturk.data = [s1.psiturk.data;s2.psiturk.data];
data.psiturk.subids = [s1.psiturk.subids;s2.psiturk.subids];
data.meadows.color.meta = [s1.meadows.color.meta;s2.meadows.color.meta];
data.meadows.color.ids = [s1.meadows.color.ids;s2.meadows.color.ids];
data.meadows.large.meta = [s1.meadows.large.meta;s2.meadows.large.meta];
data.meadows.large.ids = [s1.meadows.large.ids;s2.meadows.large.ids];

%% save data structure
cd([dr.top,dr.exp,'/Data'])
save('data_structure_full_n92.mat','data')