function [arranged_rdm, arranged_names, with_cats_rats] = arrange_rdm_full_exp2(x,filenames,item_ind,sortby)

%% setup
% if input is vector, reshape to matrix
if isvector(x)
    x = squareform(x);
end

% extract filenames
headers = {};
for i_nom = 1:length(filenames)
    headers{1,i_nom} = filenames(i_nom);
end

% add index rows/columns to RDM for sorting
x(:,end+1:end+5) = 0;
x(end+1:end+5,:) = 0;


%% add category and rating info to RDM
for i_dis = 1:length(headers)
    for i_tem = 1:size(item_ind,1)
        
        % match image name with info in item_ind
        if strcmp(cell2mat(headers{1,i_dis}),item_ind.filename{i_tem})
            
            % add category and rating info to new colums
            x(i_dis,end) = item_ind.large(i_tem);
            x(i_dis,end-1) = item_ind.color(i_tem);
            x(i_dis,end-2) = item_ind.cat2(i_tem);
            x(i_dis,end-3) = item_ind.cat1(i_tem);
            
            % and to new rows
            x(end,i_dis) = item_ind.large(i_tem);
            x(end-1,i_dis) = item_ind.color(i_tem);
            x(end-2,i_dis) = item_ind.cat2(i_tem);
            x(end-3,i_dis) = item_ind.cat1(i_tem);
            
            % add item number to headers for keeping track
            headers{2,i_dis} = item_ind.itemnr_new(i_tem);
            
        end
    end
end

% now add item numbers to RDM
x(1:end-5,end-4) = cell2mat(headers(2,:))';
x(end-4,1:end-5) = cell2mat(headers(2,:));
x(end-4:end,end-4:end) = 999;  % for sorting


%% sort
x = sortrows(x,43);
x = x';
x = sortrows(x,43);
x = x';



%% reconcile image names with now sorted RDM
for i_tag = 1:size(x,2)-5
    ind = x(i_tag,end-4);
    for i_mor = 1:size(x,2)-5
        if headers{2,i_mor} == ind
            tags{i_tag,1} = headers{1,i_mor};
        end
    end
end


%% output
arranged_rdm   = x(1:end-4,1:end-5);
arranged_names = tags;
with_cats_rats = x;


end

