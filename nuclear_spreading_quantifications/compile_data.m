% This script loads .tif files exported from Leica's LAS software, extracts the nuclear cloud and embryo mask,
% and calculates the ratio. This is the first script in the sequence of quantifying nuclear spreading.


% manually entered parameters
experiment_names = {'w1118'; 'shkl163_shkl130'; 'shkl130_EY21463'; 'pBabr_cul5_shkl163_shkl130'; 'src64KO_shkl163_shkl130'}; % string for the date of the experiment desired

compilation_ratios = cell(numel(experiment_names),1);
for i = 1:numel(experiment_names)
    j = 1;
    px_sizes = 0.727;
    while (true)

        % loading the image
        if (j < 10)
            embryo_num = strcat('00',num2str(j));
        elseif (j < 100)
            embryo_num = strcat('0',num2str(j));
        else
            embryo_num = num2str(j);
        end
        img_name = strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'.tif');
        if (exist(img_name) == 2)
            data_image = imread(img_name);
        else
            break
        end

        % make embryo mask and get the length
        if (~(exist(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_embryo_mask.mat')) == 2))
            disp('draw a polygon around the whole embryo');
            fig = figure;
            imshow(padarray(data_image*5,[100,100],0,'both'));
            h = impoly;
            embryo_mask = createMask(h);
            save(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_embryo_mask.mat'),'embryo_mask');
            close(fig);
        else
            load(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_embryo_mask.mat'));
        end
        s = regionprops(embryo_mask,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
        embryo_length = s.MajorAxisLength;

        % getting the number of nuclei
        if (~(exist(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_nuclear_pts.mat')) == 2))
            disp('click on each of the nuclei, then press enter');
            fig = figure;
            imshow(data_image);
            [pts_x,pts_y] = getpts(fig);
            pts = [pts_x,pts_y];
            save(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_nuclear_pts.mat'),'pts');
            close(fig);
        else
            load(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_nuclear_pts.mat'));
        end
        
        % getting the nuclear mask from the points
        if (~(exist(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_nuclear_mask.mat')) == 2))
            c = convhull(pts);
            fig = figure('visible','off');
            imshow(data_image);
            h = impoly(gca,pts(c,:));
            nuclear_mask = createMask(h);
            save(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_nuclear_mask.mat'),'nuclear_mask');
            close(fig);
        else
            load(strcat('exp_data/',experiment_names{i},'/','embryo_',embryo_num,'_nuclear_mask.mat'));
        end
        s = regionprops(nuclear_mask,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
        nuclei_length = s.MajorAxisLength;
        
        
        compilation_ratios{i}(j,1) = nuclei_length/embryo_length;
        compilation_ratios{i}(j,3) = nuclei_length*px_sizes;
        compilation_ratios{i}(j,4) = embryo_length*px_sizes;
        compilation_ratios{i}(j,2) = size(pts,1);
        j = j + 1;
    end
end

% transform number of nuclei to cell cycle
compilation_ratios_binned = cell(numel(compilation_ratios),1);
for i = 1:numel(compilation_ratios)
    compilation_ratios_binned{i} = NaN(size(compilation_ratios{i},1),2);
    for j = 1:size(compilation_ratios{i},1)
        compilation_ratios_binned{i}(j,2) = compilation_ratios{i}(j,1);
        compilation_ratios_binned{i}(j,1) = nextpow2(compilation_ratios{i}(j,2))+1;
        compilation_ratios_binned{i}(j,3) = compilation_ratios{i}(j,3);
        compilation_ratios_binned{i}(j,4) = compilation_ratios{i}(j,4);
    end
    compilation_ratios_binned{i} = sortrows(compilation_ratios_binned{i},1);
end

save('workspaces/compilation_ratios.mat','compilation_ratios','-v7.3');
save('workspaces/compilation_ratios_binned.mat','compilation_ratios_binned','-v7.3');