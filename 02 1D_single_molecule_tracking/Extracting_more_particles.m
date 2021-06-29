% to use this just change the indexes
%% The first part is to segment particles that were missed in the first round
for i = 1:length(actualIndex_ATPgammaS)
    if i ==1
        counter = 1;
        uniqueindexes(counter)= i;
        counter = counter+1;
    elseif length(FINAL_DATA(actualIndex_ATPgammaS(i-1)).full_length_kymo)==length(FINAL_DATA(actualIndex_ATPgammaS(i)).full_length_kymo)
    else
        uniqueindexes(counter)= i;
        counter = counter+1;
    end
end


for i = 1:length(uniqueindexes)
    if i == length(uniqueindexes)
        for j = uniqueindexes(i):length(actualIndex_ATPgammaS)
        FINAL_DATA(actualIndex_ATPgammaS(j)).kymoID = i; 
        end
    else
        for j = uniqueindexes(i):uniqueindexes(i+1)-1
        FINAL_DATA(actualIndex_ATPgammaS(j)).kymoID = i; 
        end
    end
end
 
counter = 1;
for i = 1:length(uniqueindexes)
    kymois = flip(FINAL_DATA(actualIndex_ATPgammaS(uniqueindexes(i))).full_length_kymo);
    imagesc(kymois)
    ax = gca;
    colormap(map1)
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    hold on
    if i == length(uniqueindexes)
        for j = actualIndex_ATPgammaS(uniqueindexes(i)):max(actualIndex_ATPgammaS)
        holdingplace1 = FINAL_DATA(j).crop_coords(4)-FINAL_DATA(j).crop_coords(3);
        holdingplace2 = FINAL_DATA(j).crop_coords(2)- FINAL_DATA(j).crop_coords(1);
        rectangle('Position',[FINAL_DATA(j).crop_coords(1),FINAL_DATA(j).crop_coords(3), holdingplace2 , holdingplace1], 'EdgeColor','y'); 
        end
    else
        for j = actualIndex_ATPgammaS(uniqueindexes(i)):(actualIndex_ATPgammaS(uniqueindexes(i+1))-1)
        holdingplace1 = FINAL_DATA(j).crop_coords(4)-FINAL_DATA(j).crop_coords(3);
        holdingplace2 = FINAL_DATA(j).crop_coords(2)- FINAL_DATA(j).crop_coords(1);
        rectangle('Position',[FINAL_DATA(j).crop_coords(1),FINAL_DATA(j).crop_coords(3), holdingplace2 , holdingplace1], 'EdgeColor','y'); 
        end
    end       
    hold off
    ax.Color = 'none';
    set(gca,'Visible','off')
%     ax.DataAspectRatio = [1 1 1];
    caxis([0 8])
    
    ksize = size(kymois);
    kw=ksize(:,2);
    kh=ksize(:,1);
    
    
    prompt = 'How many traces to extract?';
    num_of_traces = input(prompt);
    for j = 1:num_of_traces
        p = ginput(2); 
        sp(1) = min(floor(p(1)), floor(p(2))); %xmin
        sp(2) = min(floor(p(3)), floor(p(4))); %ymin
        sp(3) = max(ceil(p(1)), ceil(p(2)));   %xmax
        sp(4) = max(ceil(p(3)), ceil(p(4)));   %ymax
        if sp(1)<1
            sp(1)=1;
        end
        if sp(2)<1
            sp(2)=1;
        end
        if sp(3)>(kw-1)
            sp(3)=kw;
        end
        if sp(4)>(kh-1)
            sp(4)=kh;
        end
        crop = kymois(sp(2):sp(4), sp(1): sp(3),:);
        structure_name(counter).name = FINAL_DATA(actualIndex_ATPgammaS(uniqueindexes(i))).mol_id;
        structure_name(counter).crop = crop;
        structure_name(counter).crop_coordinates = [sp(2),sp(4),sp(1),sp(3)];
        crop_fig = figure;
          if max(kymois,[], 'all') > 30
            imagesc(crop,[0,max(crop, [], 'all')/2])  
            set(crop_fig,'Position',[594.0000, 33.3333, 334.0000, 612.6667]);
            disp(['max intensity is ', num2str(max(kymois,[], 'all'))]);
          else
            imagesc(crop);      
            set(gcf,'Position',[594.0000, 33.3333, 334.0000, 612.6667]);
          end
        ksize = size(crop);
        kw=ksize(:,2);
        kh=ksize(:,1);
        p = ginput(2); 
        sp(1) = min(floor(p(1)), floor(p(2))); %xmin
        sp(2) = min(floor(p(3)), floor(p(4))); %ymin
        sp(3) = max(ceil(p(1)), ceil(p(2)));   %xmax
        sp(4) = max(ceil(p(3)), ceil(p(4)));   %ymax
            if sp(1)<1
                sp(1)=1;
            end
            if sp(2)<1
                sp(2)=1;
            end
            if sp(3)>(kw-1)
                sp(3)=kw;
            end
            if sp(4)>(kh-1)
                sp(4)=kh;
            end
        BG = crop(sp(2):sp(4), sp(1): sp(3),:);
        BG = reshape(BG.',1,[]);
        d=mean(BG);
        dsd=sqrt(2)*std(BG);%% fix a mistake in the or
        dmin=d-(2*dsd);
        dmax=d+(2*dsd);
        structure_name(counter).dmin = dmin;
        structure_name(counter).dmax = dmax;        
        structure_name(counter).d = d;
        crop_fig = figure;
        imagesc(crop,[0,max(crop, [], 'all')/2]);
        set(crop_fig,'Position',[594.0000, 33.3333, 334.0000, 612.6667]); 
        p = ginput(2); 
        sp(11) = min(floor(p(1)), floor(p(2))); %xmin
        sp(12) = min(floor(p(3)), floor(p(4))); %ymin
        sp(13) = max(ceil(p(1)), ceil(p(2)));   %xmax
        sp(14) = max(ceil(p(3)), ceil(p(4)));   %ymax
            if sp(11)<1
                sp(11)=1;
            end
            if sp(12)<1
                sp(12)=1;
            end
            if sp(13)>(kw-1)
                sp(13)=kw;
            end
            if sp(14)>(kh-1)
                sp(14)=kh;
            end
        peak = crop(sp(12):sp(14), sp(11): sp(13),:);
        peak = reshape(peak.',1,[]);
        a=mean(peak);
        amin=a-(3*dsd);
        amax=a+(3*dsd);
        structure_name(counter).amin = amin;
        structure_name(counter).amax = amax;
        structure_name(counter).a = a;
        sz=size(crop);
        t=sz(:,1);
        w=sz(:,2);
        b=(sp(11)+sp(13))/2;
        bmin=b-5;
        bmax=b+5;
        structure_name(counter).bmin = bmin;
        structure_name(counter).bmax = bmax;
        structure_name(counter).b = b;
        wrange=[1:1:w];
        transpose(wrange);
        structure_name(counter).wrange = wrange;
        clf
        counter = counter + 1;
        
        kymofig = figure;
        if j ~= num_of_traces
            kymois = flip(FINAL_DATA(actualIndex_ATPgammaS(uniqueindexes(i))).full_length_kymo);
            imagesc(kymois)
            ax = gca;
            colormap(map1)
            ax.XTickLabel = [];
            ax.YTickLabel = [];
            hold on
            if i == length(uniqueindexes)
                for j = actualIndex_ATPgammaS(uniqueindexes(i)):max(actualIndex_ATPgammaS)
                holdingplace1 = FINAL_DATA(j).crop_coords(4)-FINAL_DATA(j).crop_coords(3);
                holdingplace2 = FINAL_DATA(j).crop_coords(2)- FINAL_DATA(j).crop_coords(1);
                rectangle('Position',[FINAL_DATA(j).crop_coords(1),FINAL_DATA(j).crop_coords(3), holdingplace2 , holdingplace1], 'EdgeColor','y'); 
                end
            else
                for j = actualIndex_ATPgammaS(uniqueindexes(i)):(actualIndex_ATPgammaS(uniqueindexes(i+1))-1)
                holdingplace1 = FINAL_DATA(j).crop_coords(4)-FINAL_DATA(j).crop_coords(3);
                holdingplace2 = FINAL_DATA(j).crop_coords(2)- FINAL_DATA(j).crop_coords(1);
                rectangle('Position',[FINAL_DATA(j).crop_coords(1),FINAL_DATA(j).crop_coords(3), holdingplace2 , holdingplace1], 'EdgeColor','y'); 
                end
            end       
            hold off
            ax.Color = 'none';
            set(gca,'Visible','off')
        %     ax.DataAspectRatio = [1 1 1];
            caxis([0 8])

            ksize = size(kymois);
            kw=ksize(:,2);
            kh=ksize(:,1);
        else
        end
    end
    close all
end

%% The second part is to fit to determine the position of the particle

for i = 1:length(structure_name_gammaS)
%     structure_name_gammaS(i).crop = rot90(structure_name_gammaS(i).crop);
    sz=size(structure_name_gammaS(i).crop);
    t=sz(:,1);
    tt=sz(:,2);
    gaussfitting(i).name = structure_name_gammaS(i).name;
    fittingcell = {};
    gof_results = {};
    datacell = {};
    height = [];
    width = [];
    pixel = [];
    timepix = [];
    intensity = [];
    for j=[1:t]
            subcrop=structure_name_gammaS(i).crop(j,:);
            transpose(subcrop);  
            fitresult = cell( 2, 1 );
            gof = struct( 'sse', cell( 2, 1 ), ...
            'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
            [xData, yData] = prepareCurveData([1:1:tt], subcrop);
            ft = fittype( 'd + (a*exp(-((x-b)/sqrt(2)*c)^2))', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [structure_name_gammaS(i).amin structure_name_gammaS(i).bmin -2 structure_name_gammaS(i).dmin];
            opts.StartPoint = [structure_name_gammaS(i).a structure_name_gammaS(i).b -1.3 structure_name_gammaS(i).d];
            opts.Upper = [structure_name_gammaS(i).amax structure_name_gammaS(i).bmax -1 structure_name_gammaS(i).dmax];
            [fitresult, gof] = fit( xData, yData, ft, opts );
            fitcoef=coeffvalues(fitresult);
            height(j)=fitcoef(:,1);
            width(j)=abs(2*fitcoef(:,3));
            pixel(j)=fitcoef(:,2);
            timepix(j)=j;
            intensity(j)= height(j)*width(j); %fitcoef(:,1)*2*(fitcoef(:,3);
            structure_name_gammaS(i).a=fitcoef(:,1);
            dsd=(structure_name_gammaS(i).d - structure_name_gammaS(i).dmin)/2;
            structure_name_gammaS(i).amin=structure_name_gammaS(i).a-dsd;
            structure_name_gammaS(i).amax=structure_name_gammaS(i).a+dsd;
            structure_name_gammaS(i).b=fitcoef(:,2);
            structure_name_gammaS(i).bmin=structure_name_gammaS(i).b-5;
            structure_name_gammaS(i).bmax=structure_name_gammaS(i).b+5;
            fittingcell{j} = fitresult;
            gof_results{j} = gof;
            datacell{j} = [xData, yData];  
    end
    disp([num2str(j), '/', num2str(t)])
    save_these_results.fitting = fittingcell;
    save_these_results.gof = gof_results;
    save_these_results.data = datacell;
    save_these_results.height = height;
    save_these_results.width = width;
    save_these_results.pixel = pixel;
    save_these_results.timepix = timepix;
    save_these_results.intensity = intensity;
    gaussfitting(i).data = save_these_results; 
end

%%
i = 8;
crop = rot90(structure_name_gammaS(i).crop); 
imagesc(crop)
hold on
plot(gaussfitting(i).data.pixel,gaussfitting(i).data.timepix,'y')
hold off
