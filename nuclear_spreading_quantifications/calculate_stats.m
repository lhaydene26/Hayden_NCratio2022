% This script loads data generated with the compile_data.m script, transforms the structure of the data,
% and performs the K-S statistical test for signficiance. This is the second script in the sequence of 
% quantifying nuclear spreading.

% loads the data
load('workspaces/compilation_ratios_binned.mat');

% transform binned cell cycle versus ratio data into name x cc x ratio data
compilation_ratios_cc = cell(numel(compilation_ratios_binned),5);
for i = 1:numel(compilation_ratios_binned)
    for j = 1:5
        l = 1;
        for k = 1:size(compilation_ratios_binned{i},1)
            if (compilation_ratios_binned{i}(k,1) == j+3)
                compilation_ratios_cc{i,j}(l,1) = compilation_ratios_binned{i}(k,2);
                compilation_ratios_cc{i,j}(l,2) = compilation_ratios_binned{i}(k,3);
                compilation_ratios_cc{i,j}(l,3) = compilation_ratios_binned{i}(k,4);
                l = l + 1;
            end
        end
        if (numel(compilation_ratios_cc{i,j}) > 0)
            clear compilation_ratios_cc2
            compilation_ratios_cc2(:,1) = compilation_ratios_cc{i,j}(~isoutlier(compilation_ratios_cc{i,j}(:,1)),1);
            compilation_ratios_cc2(:,2) = compilation_ratios_cc{i,j}(~isoutlier(compilation_ratios_cc{i,j}(:,1)),2);
            compilation_ratios_cc2(:,3) = compilation_ratios_cc{i,j}(~isoutlier(compilation_ratios_cc{i,j}(:,1)),3);
            compilation_ratios_cc{i,j} = compilation_ratios_cc2;
        end
    end
end
save('workspaces/compilation_ratios_cc','compilation_ratios_cc','-v7.3');

% performs statistical test
% rows are conditions, columns are 1) ratio, 2) nuclear cloud length, and 3) embryo length
cc_num = 6;
wt_stat_ks = NaN(size(compilation_ratios_cc,1),3);
shkl_stat_ks = NaN(size(compilation_ratios_cc,1),3);
for i = 1:size(compilation_ratios_cc,1)
    for j = 1:3
        if (numel(compilation_ratios_cc{i,cc_num-3}(:,j)) > 4)
            wt_stat_ks(i,j) = kstest2(compilation_ratios_cc{i,cc_num-3}(:,j),compilation_ratios_cc{1,cc_num-3}(:,j));
            shkl_stat_ks(i,j) = kstest2(compilation_ratios_cc{i,cc_num-3}(:,j),compilation_ratios_cc{2,cc_num-3}(:,j));
        end
    end
end
