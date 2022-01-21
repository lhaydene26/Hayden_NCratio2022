function sensorActivity = funct_plot_RhoActivity(data)
    % This function accepts a data structure containing file location, experiment name,
    % and the specific frames to use from the movie, averages the data over the embryo cortex,
    % and generates a plot of the sensor data over time

    sensorActivity = cell(numel(data.loc),1);
    for i = 1:numel(data.loc)
        fname = strcat(data.loc{i}, data.exp_name{i}, '_ch00_stack.tif');

        info = imfinfo(fname); %gets infos on stack
        tmax = numel(info); %gets number of images in stack
        
        im = cell(tmax,1);
        for j = 1:tmax
            im{j} = im2double(imread(fname, j, 'Info', info)); %reads in image
        end

        % embryo mask
        if (exist(strcat(data.loc{i}, 'mask_embryo.mat'),'file') == 2)
            load(strcat(data.loc{i}, 'mask_embryo.mat'));
        else
            fig = figure(1);
            imagesc(im{1})
            axis equal
            h = impoly;
            mask_embryo = createMask(h);
            save(strcat(data.loc{i}, 'mask_embryo.mat'),'mask_embryo','-v7.3')
            close(fig);
        end  
        
        sensorActivity{i} = NaN(tmax-1,1);
        for j = 1:tmax-1
            sensorActivity{i}(j) = nanmean(im{j}(mask_embryo));
        end
        sensorActivity{i} = sensorActivity{i}(data.frames{i});
        sensorActivity{i} = sensorActivity{i}/min(sensorActivity{i});
        sensorActivity{i} = sgolayfilt(sensorActivity{i},3,15);
    end

    fig = figure;
    hold on;
    for i = 1:numel(sensorActivity)
        plot(sensorActivity{i})
    end
    
end