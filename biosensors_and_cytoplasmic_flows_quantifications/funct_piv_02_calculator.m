function funct_piv_02_calculator(data_loc, exp_name)
    % This function accepts a file location and name and performs PIV

    % optional checks
    check_filter=0;
    check_points=0;
    fig=0;

    fname = strcat(data_loc, exp_name, '_cyto_mean.tif');

    info = imfinfo(fname); %gets infos on stack
    tmax = numel(info); %gets number of images in stack
    im{1} = imread(fname, 1, 'Info', info); %reads in image z-th
    im{1}=fliplr(im{1});
    ss = size(im{1});

    h = fspecial('gaussian',5,1); %% changing this value changes the thresholds in correlation
    for i = 1:tmax 
        im{i} = imread(fname, i, 'Info', info); %reads in image z-th
        im{i}=fliplr(im{i});
        if check_filter
            figure(1)
            imagesc(im{i})
        end
        im{i} =  imfilter(im{i},h);
        if check_filter
            figure(2)
            imagesc(im{i})
            return
        end
        if check_points
            figure(1)
            hold all
            imagesc(im{i});
        end
    end

    deltax = 35; %% 35 better than 30
% deltax = 30;
    deltay = 35;
    ws = 35+deltax;
    hs = 20+deltay; % how much images move depend on dt as well dx = v * dt; hs should be smthing like 3*dx
% hs = 27+deltay;

    %identifies background
    r = im{1};
    ixbg=find(r>500);
    scra  = zeros(size(r));
    scra(ixbg)=1;
    scra = bwareaopen(scra,20);
    scra = imclose(scra,ones(2));
    scra = imfill(scra,'holes');
    scra(1,:) = 0;
    scra(end,:) = 0;
    scra = imerode(scra,ones(40));
    if (sum(scra(:)) == 0)
        error('background not detected');
    end

    ixbg = find(~scra);
    [ybg,xbg]=ind2sub(size(r),ixbg);

    if check_points
        figure(3)
        hold all
        plot(xbg,ybg,'bs')
        xlim([1 ss(2)])
        ylim([1 ss(1)])
    end
    %generate random points
    npmax=1000;
    xr=[];yr=[];
    np=numel(xr);
    while np<npmax
        xr = [xr; ws+floor((ss(2)-2*ws)*rand(npmax,1))]; %% increased the region where points are taken
        yr = [yr; hs+floor((ss(1)-2*hs)*rand(npmax,1))];
        [xx]=setdiff([xr yr],[xbg,ybg],'rows');
        xx = unique(xx,'rows');
        xr = xx(:,1);
        yr = xx(:,2);
        np=numel(xr);
    end
    if check_points
        figure(3)
        plot(xr,yr,'rs')
        figure(1)
        hold all
        plot(xr,yr,'rs')
        return
    end
%     np = numel(xr)

    offx = cell(1,tmax);
    offy = cell(1,tmax);
    xt = cell(1,np);
    yt = cell(1,np);

    ddx = 1; %change if frames are not temporally uniform
    dt=6; %was 3
    for t = dt+1:tmax  % use all frames
        q=0;
        tstart=t-dt;

        xs=[];ys=[];
        dxm = [];dym = [];
        if fig
            figure(2)
            hold all
            imagesc(im{t-dt}), colormap(gray(255))
            hold all
            axis equal

        end
        for i = 1:np
            cc = [];
            xc = xr(i);
            yc = yr(i);

            imb = im{t-dt};
            ima = im{t}; %normal
            recta=[ xc-round(ws/2) yc-round(hs/2) ws hs];
            ia=imcrop(ima,recta);
            rectb = [xc-round(deltax/2) yc-round(deltay/2) deltax deltay];
            ib=imcrop(imb,rectb);
            if std(im2double(ib(:))) == 0
                display('no std');
                continue
            end
            c=normxcorr2(ib,ia);
            % offset found by correlation
            [max_c, imax] = max(c(:));

            if fig
                figure(10)
                imagesc(ia)
                axis equal
                figure(11)
                imagesc(ib)
                axis equal
                figure(12)
                surf(c)
                max_c
                pause
            end

            if max_c >0.5
                [ypeak, xpeak] = ind2sub(size(c),imax(1));
                corr_offset = [(xpeak-size(ib,2)) (ypeak-size(ib,1))];

                % relative offset of position of subimages
                rect_offset = [(recta(1)-rectb(1)) (recta(2)-rectb(2))];
                % total offset
                offset = corr_offset + rect_offset;
                if (offset(1)^2+offset(2)^2<0.1*ws^2) %max displacements should not exceed a threshold.
                    x = [xc xc+offset(1)];
                    y = [yc yc+offset(2)];
                    q = q+1;
                    dxm(q) = (x(2)-x(1))/ddx;
                    dym(q) = (y(2)-y(1))/ddx;

                    if isnan(dxm(q))
                        q
                        return
                    end
                    if isnan(dym(q))
                        q
                        return
                    end

                    xs(q) = x(1); %%% these last 7 lines should be here
                    ys(q) = y(1);

                else
                    continue
                end
            else
                continue
            end
            %             waitbar(q/np,h,'On my way...')


        end
%         numel(dxm)

        scale=1/dt;
        dxm=scale*dxm;
        dym=scale*dym;

        s=size(dxm);

%         dxm_all(t,:)=dxm;


        dxm_all{t,1}=dxm;
        dym_all{t,1}=dym;


        dxm_plot=dxm; dxm_plot(isnan(dxm))=0;
        dym_plot=dym; dym_plot(isnan(dym))=0;


        ix = find(dxm_plot==0 & dym_plot==0);


        if fig
            figure(2)
            plot_velocity=quiver(xs,ys,dxm_plot,dym_plot,0,'color','r','linewidth',0.35);
            axis equal
            xlim([1 ss(2)])
            ylim([1 ss(1)])

            hold on
            plot(xs(ix),ys(ix),'y.','markersize',1)


            %h=figure(2)
            figname = sprintf('%svelocity_field_%03i.tif',fname,t-dt);
            axis off
            %set(gcf, 'PaperPositionmode', 'auto');
            %print('-dtiff',figname,'-r200')
            clf(2);
        end

        xs_all{t,1}=xs;
        ys_all{t,1}=ys;

        figname = strcat(fname(1:end-23), 'velocity_fields2\',sprintf('velocity_field_%03i.mat',t-dt));
        save(figname,'xs','ys','dxm','dym');
        t/tmax
    end
end