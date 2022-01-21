function data_NCratio_EDprob = funct_NCratio_EDprob(dirs, gs, dx, nboxes_AP, nboxes_DV)
    % This function accepts the 1) file locations for embryos to process, 2) community radii
    % in which to calculate the N/C ratio, 3) pixel size for movies (standard: 0.56), 4) number of grid points
    % across the AP axis (standard: 80), and 5) number of grid points across the DV axis (standard: 30). 
    %
    % This function processes each embryo and calculates the N/C ratio in circles of different community
    % radii
    
    [X,Y] = meshgrid(linspace(1, 800, nboxes_AP), linspace(1, 300, nboxes_DV));

    % for each of the conditions
    data_NCratio_EDprob = NaN(numel(gs), 1, 2);
    storage_index = ones(numel(gs),1);
    for i = 1:numel(dirs)
        % for each embryo in the folder
        Folder = dirs{i};
        filePattern = fullfile(Folder, '*.tif');
        tifFiles   = dir(filePattern);
        k = 1;
        for j = 1:length(tifFiles)
            temp = tifFiles(k);
            if (temp.name(end-4) == 's')
                tifFiles(k) = [];
            else
                k = k + 1;
            end
        end
        for j = 1:length(tifFiles)
            % getting file name
            baseFileName = tifFiles(j).name;
            fullFileName = fullfile(Folder, baseFileName);

            % load the embryo data
            data = imread(fullFileName);
            
            % load the embryo and region masks
            mask_embryo = logical(imread(strcat(dirs{i}, 'masks\', tifFiles(j).name(1:end-4),'_mask.tif')));
            if (exist(strcat(dirs{i}, 'masks\', tifFiles(j).name(1:end-4),'_ND.tif'),'file') == 2)
                mask_ND = logical(imread(strcat(dirs{i}, 'masks\', tifFiles(j).name(1:end-4),'_ND.tif')));
                mask_ED = and(mask_embryo,~mask_ND);
            else
                mask_ND = NaN;
                mask_ED = NaN;
            end  
        
            % segment the nuclei
            data_segmented = funct_segment_nuclei(data, mask_embryo, 0.35);
            
            % extract nuclei positions from segmented data
            stats = regionprops(data_segmented, 'centroid');
            nuclei_pts = NaN(numel(stats),2);
            temp = [stats.Centroid];
            m = 1;
            for k = 1:2:(numel(temp)-1)
                nuclei_pts(m,1) = temp(k);
                nuclei_pts(m,2) = temp(k+1);
                m = m + 1;
            end
            
            % save nuclei pts
            save(strcat(fullFileName(1:end-4),'_nuclei_pts.mat'),'nuclei_pts','-v7.3');
            
            % for each size box
            for m = 1:numel(gs)
                for n = 1:(nboxes_AP*nboxes_DV)
                    if (~mask_embryo(round(Y(n)),round(X(n))))
                        continue;
                    else
                        N = rangesearch(nuclei_pts, [X(n), Y(n)], gs(m)/dx);
                        N = numel(N{1});

                        [XX, YY] = meshgrid(1:size(data,2), 1:size(data,1)); % grid for the x and y coords
                        dRegion = gs(m)/dx; % radius
                        boxRegion = sqrt((XX-X(n)).^2 + (YY-Y(n)).^2) < dRegion;
                        boxRegion = and(boxRegion,mask_embryo);
                        C = sum(boxRegion,'all')*dx^2;
  
                        data_NCratio_EDprob(m,storage_index(m),1) = N/C;
                        
                        data_NCratio_EDprob(m,storage_index(m),3) = X(n);
                        data_NCratio_EDprob(m,storage_index(m),4) = Y(n);

                        if (strcmp(dirs{i}(15:18), 'cc13'))
                            data_NCratio_EDprob(m, storage_index(m), 2) = 1;
                        elseif (strcmp(dirs{i}(15:18), 'cc14'))
                            data_NCratio_EDprob(m, storage_index(m), 2) = 0;
                        elseif (strcmp(dirs{i}(15:19), 'split'))
                            if (mask_ND(round(Y(n)),round(X(n))))
                                data_NCratio_EDprob(m, storage_index(m), 2) = 0;
                            elseif (mask_ED(round(Y(n)),round(X(n))))
                                data_NCratio_EDprob(m, storage_index(m), 2) = 1;
                            else
                                continue;
                            end                         
                        else
                            error('folder naming error');
                        end
                        storage_index(m) = storage_index(m) + 1;
                    end
                end
            end
            i
            j
        end
    end
end