function data_cc14_NCratio_avg = funct_calc_NCratio_cc14(data_cc14_loc)
    % This function accepts the embryo data which halt normally at cc14 (specifically the embryo
    % identification numbers) and calculates the average NC ratio across the embryo (to be used to normalize future N/C
    % ratio values
    
    load(data_cc14_loc);
    data_cc14_NCratio_avg = NaN(size(data_NCratio_EDprob,1),1);
    for i = 1:size(data_NCratio_EDprob,1)
        data_cc14_NCratio_avg(i) = nanmean(data_NCratio_EDprob(i,:,1));
    end
end