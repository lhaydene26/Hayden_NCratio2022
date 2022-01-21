function [data_averaged, data_fit_params] = funct_binANDavg(data_NCratio_EDprob, data_cc14_NCratio_avg, NCratio_bins)
    % This function accepts the processed N/C ratio data, average N/C ratio data, and the edges
    % of N/C ratio bins to use (21, from 0 to 0.04). The script bins to data and then calculates
    % the average in each bin
    
    % bin the data for averaging
    data_binned = cell(size(data_NCratio_EDprob,1),numel(NCratio_bins));
    for i = 1:size(data_NCratio_EDprob,1)
        for j = 1:size(data_NCratio_EDprob,2)
            y = discretize(data_NCratio_EDprob(i,j,1), NCratio_bins);
            if (~isnan(y))
                data_binned{i,y}(end+1,1) = data_NCratio_EDprob(i,j,1);
                data_binned{i,y}(end,2) = data_NCratio_EDprob(i,j,2);
            end
        end
    end

    data_averaged = NaN(size(data_NCratio_EDprob,1),numel(NCratio_bins),3);
    for i = 1:size(data_NCratio_EDprob,1)
        for j = 1:numel(NCratio_bins)
            if (numel(data_binned{i,j}) > 0)
                data_averaged(i,j,1) = nanmean(data_binned{i,j}(:,1));
                if (~isnan(data_averaged(i,j,1)))
                    data_averaged(i,j,2) = nanmean(data_binned{i,j}(:,2));
                    data_averaged(i,j,3) = nanstd(data_binned{i,j}(:,2));
                end
            end
        end
    end
    
    
    % fit a logistic curve to each point and plot the 
    steepness = NaN(size(data_averaged,1),1);
    pos = steepness;
    ft = fittype('1/(1+exp(-(-x-b)/k))');
    for i = 1:size(data_averaged,1)
    %     if (i>1)
    %         f = fit(data_averaged(i,~isnan(data_averaged(i,:,2)),1)',data_averaged(i,~isnan(data_averaged(i,:,2)),2)',ft,'StartPoint', [f.b, f.k]);
    %     else
            f = fit((data_averaged(i,~isnan(data_averaged(i,:,2)),1)/data_cc14_NCratio_avg(i))',data_averaged(i,~isnan(data_averaged(i,:,2)),2)',ft,'StartPoint', [0.015, 400]);
    %     end
        steepness(i) = f.k;
        pos(i) = f.b;
    end
    
    data_fit_params.pos = pos;
    data_fit_params.steepness = steepness;
    
end
