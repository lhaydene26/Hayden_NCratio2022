function data_segmented = funct_segment_nuclei(data, mask_embryo, p)
    % This is an accessory function called by funct_NCratio_EDprob.m to segment nuclei
    
    b = imgaussfilt(data,2);
    T = adaptthresh(b,p,'ForegroundPolarity','bright');
    data_segmented = imbinarize(b,T);
    for i = 1:size(data_segmented,1)
        for j = 1:size(data_segmented,2)
            if (~mask_embryo(i,j))
                data_segmented(i,j) = 0;
            end
        end
    end

    
    
end