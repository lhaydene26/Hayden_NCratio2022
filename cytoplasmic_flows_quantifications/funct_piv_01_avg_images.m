function funct_piv_01_avg_images(data_loc, exp_name)
    % This function accepts a file location and name and smoothes the cytoplasmic data through averaging for
    % use in PIV
    fname = strcat(data_loc, exp_name, '_cyto.tif');

    info = imfinfo(fname); %gets infos on stack
    tmax = numel(info); %gets number of images in stack

    for i = 1:tmax
        im(:,:,i) = im2double(imread(fname, i, 'Info', info)); %reads in image z-th
    end
    q=0;
    newfilename = strcat(fname(1:end-4), '_mean.tif');
    di=1;
    for i = 1+di:tmax-1-di
        q = q+1;
        newim{q} = im2uint16(mean([im(:,:,i-di:i+di)],3));
        if q==1
            imwrite(newim{q}, newfilename,'compression','none');
        else
            imwrite(newim{q}, newfilename, 'WriteMode', 'append','compression','none');
        end
        i
    end

end