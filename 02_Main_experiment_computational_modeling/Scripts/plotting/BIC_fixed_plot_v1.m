function BIC_fixed_plot_v1( bicF, title_tag )
order = [1:5];
nn = size( bicF, 2);
bicF_rev = zeros(1,nn);


for ii = 1:nn; 
    bicF_rev(1,ii) = bicF(1,nn+1-ii);
end

for ii = 1:nn; 
    bicF_ord(1,ii) = bicF(1,order(ii));
end
for ii = 1:nn;
bicF_rev_ord(1,ii) = bicF_ord(1,nn+1-ii);
end
bicF_rev_ord = bicF_rev_ord - bicF_rev_ord(4);


%bicF_rev = bicF_rev - bicF_rev(5);
figure
barh( bicF_rev_ord, 'k' )
%xlim([-380,400]);
title( title_tag )

