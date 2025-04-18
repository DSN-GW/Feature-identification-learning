function [arranged_rdm,arranged_names] = arrange_rdm(x,filenames,item_ind)

%  if input is vector, reshape to matrix
if isvector(x)
    x = squareform(x);
end


%% add category info to rdm
for i_dis = 1:length(filenames)
%this here tells you the location of the pic in your rdm
 check(i_dis) = find(strcmp(filenames,item_ind.filename{i_dis}));
end 
%% add rdm to item_ind // one master file with all the info
item_ind.rdm = check';

%% reshape x mat by the correct sort 

item_ind=sortrows(item_ind,6);
y = x;
%% this is not at all elegant and we should change
y(:,121)=item_ind.itemnr_new;
y=sortrows(y,121);
y=y(:,1:120);
y(121,:)=item_ind.itemnr_new';
y=y';
y=sortrows(y,121);
y=y(:,1:120); % this should be sorted correctly


%  output
arranged_rdm = y;
arranged_names = item_ind;
%with_cats = x;


end

