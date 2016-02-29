% code for analyzing the contrast responses of a divisively normalized
% contrast encoder
clear all; close all; addpath(genpath('./helper_functions'));
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontSize', 20)

%% set up filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rf_size         = -80:0.2:80;                       % positions at which to define the filter
[x,y]           = meshgrid(rf_size,rf_size);        % x and y coordinates of each position
csig            = 3;                                % sigma of central Gaussian
center          = make_2D_gaussian(rf_size,csig);   % central Gaussian

%% make basic stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conts           = 0:0.01:1;     % Weber contrast of spot
diam            = 4*csig;       % Diameter of stimulus
c_mask          = double(((x.^2+y.^2)<=(diam/2)^2)); % make a positive value mask for stimulus

%% define scales for family of filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssigs       = [1.25 1.5 2 3 4 6];       % scale factors for surround (x center)
assigs      = [1 1.25 1.5 2 3 4 6];     % scale factors for adaptation (x center)

%% make colormaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cex     = hslcolormap(1500,0.1,0.7,[1 0.5]);
cin     = hslcolormap(abs(round(1500*(-0.0992))),0.5,0.6,[1 0.5]);
csrf    = [flipud(cin); cex];
cada    = hslcolormap(0.3,0.5,[1 0.5]);

%% plot basic filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
    
    ssig            = 2*csig;                           % sigma of antagonistic surround Gaussian
    surround        = make_2D_gaussian(rf_size,ssig);   % surround Gaussian
    srf             = center - surround;                % DOG
    ada             = surround;                         % adaptive pool
    cropamt         = 310;                              % crop edges by this much for visuals
    
    f = figure; hold on; setupfig(6,3,12);
    surf(crop_edges(srf,cropamt),'edgecolor','none');
    colormap(csrf); view([-14 22]); set(gca,'Xtick',[],'YTick',[],'ZTick',0,'Zgrid','on');
    axis([0 length(rf_size)-cropamt*2 0 length(rf_size)-cropamt*2]); hcb = colorbar; set(hcb,'YTick',[0])
    export_fig('./figures/basicfilter_1.png','-r300'); close(f);
    
    f = figure; hold on; setupfig(6,3,12);
    surf(crop_edges(ada,cropamt),'edgecolor','none');
    colormap(cada); view([-14 22]); set(gca,'Xtick',[],'YTick',[],'ZTick',[]);
    axis([0 length(rf_size)-cropamt*2 0 length(rf_size)-cropamt*2]); hcb = colorbar; set(hcb,'YTick',[0])
    export_fig('./figures/basicfilter_2.png','-r300'); close(f);
    
    f = figure; hold on; setupfig(6,6,12);
    surf(crop_edges(srf,cropamt),'edgecolor','none');
    colormap(csrf); set(gca,'Xtick',[],'YTick',[],'ZTick',0,'Zgrid','on');
    axis([0 length(rf_size)-cropamt*2 0 length(rf_size)-cropamt*2]); axis square;
    export_fig('./figures/basicfilter_3.png','-r300'); close(f);
    
    f = figure; hold on; setupfig(6,6,12);
    surf(crop_edges(ada,cropamt),'edgecolor','none');
    colormap(cada); grid off; box on; set(gca,'Xtick',[],'YTick',[],'ZTick',[]);
    axis([0 length(rf_size)-cropamt*2 0 length(rf_size)-cropamt*2]); axis square;
    export_fig('./figures/basicfilter_4.png','-r300'); close(f);
    
end

%% CRF w and w/o norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
    
    gray    = 50;       % intensity of middle gray
    gray2   = 100;      % intensity of middle gray for brighter stimuli
    stimim  = [];       % initialize matrix for stimulus images
    stimim2 = [];       % 2x brighter stimulus
    
    for k = 1:length(conts)                             % for each contrast level
        
        % regular mean luminance
        sON         = (gray.*conts(k).*c_mask) + gray;      % generate the Bright and Dark stimuli
        sOFF        = (gray.*-conts(k).*c_mask) + gray;
        
        rONn(k)     = model_response('ON',sON,srf,ada,0);   % contrast responses without norm
        rOFFn(k)    = model_response('OFF',sOFF,srf,ada,0);
        rON(k)      = model_response('ON',sON,srf,ada,1);   % contrast responses with norm
        rOFF(k)     = model_response('OFF',sOFF,srf,ada,1);
        
        % higher mean luminance
        sONg2       = (gray2.*conts(k).*c_mask) + gray2;    % generate the Bright and Dark stimuli
        sOFFg2      = (gray2.*-conts(k).*c_mask) + gray2;
        
        rONng2(k)   = model_response('ON',sONg2,srf,ada,0); % contrast responses without norm
        rOFFng2(k)  = model_response('OFF',sOFFg2,srf,ada,0);
        rONg2(k)    = model_response('ON',sONg2,srf,ada,1); % contrast responses with norm
        rOFFg2(k)   = model_response('OFF',sOFFg2,srf,ada,1);
        
        % store stimulus image
        if ~mod(conts(k),0.1); stimim = cat(2,stimim,cat(1,crop_edges(sOFF,350),crop_edges(sON,350))); end
        if ~mod(conts(k),0.1); stimim2 = cat(2,stimim2,cat(1,crop_edges(sOFFg2,350),crop_edges(sONg2,350))); end
    end
    
    % save stimulus image
    figure; hold on; setupfig(12,4,20); axis image; colormap gray; set(gca,'ytick',[],'xtick',[]);
    imagesc(stimim,[0 200]);
    export_fig('./figures/basicCRFstim.png');
    
    figure; hold on; setupfig(12,4,20); axis image; colormap gray; set(gca,'ytick',[],'xtick',[]);
    imagesc(stimim2,[0 200]);
    export_fig('./figures/basicCRFstim_brighter.png');
    
    % CRF w/o norm
    figure; hold on; setupfig(6,6,20); title('CRF without norm'); axis square;
    h(1) = plot(100*conts, rONn,'color',ColorIt('r'));
    h(2) = plot(100*conts, rOFFn,'--','color',ColorIt('b'));
    legend(h,'ON','OFF'); set(gca,'ytick',[]);
    ylim([0 max([rONng2 rOFFng2])]); xlabel('Weber contrast'); ylabel('filter response');
    export_fig('./figures/basicCRFnonorm.pdf');
    
    figure; hold on; setupfig(6,6,20); title('CRF without norm'); axis square;
    h(1) = plot(100*conts, rONng2,'color',ColorIt('r'));
    h(2) = plot(100*conts, rOFFng2,'--','color',ColorIt('b'));
    legend(h,'increments','decrements'); set(gca,'ytick',[]);
    ylim([0 max([rONng2 rOFFng2])]); xlabel('Weber contrast'); ylabel('filter response');
    export_fig('./figures/basicCRFnonorm_brighter.pdf');
    
    % CRF w/ norm
    figure; hold on; setupfig(6,6,20); title('CRF with norm'); axis square;
    plot(100*conts, rON,'color',ColorIt('r'));
    plot(100*conts, rOFF,'--','color',ColorIt('b'));
    ylim([0 max([rON rOFF])]); xlabel('Weber contrast'); ylabel('filter response'); set(gca,'ytick',[]);
    export_fig('./figures/basicCRF.pdf');

    figure; hold on; setupfig(6,6,20); title('CRF with norm'); axis square;
    plot(100*conts, rONg2,'color',ColorIt('r'));
    plot(100*conts, rOFFg2,'--','color',ColorIt('b'));
    ylim([0 max([rON rOFF])]); xlabel('Weber contrast'); ylabel('filter response'); set(gca,'ytick',[]);
    export_fig('./figures/basicCRF_brighter.pdf');
end


%% nonlinearity indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
    
    % intensity of middle gray
    gray = 50;                   
    
    % plot stimulus size
    cf = figure; hold on; setupfig(6,6,12); axis square; surf(1-c_mask,'edgecolor','none');
    colormap gray;
    set(gca,'Xtick',[],'YTick',[],'ZTick',0); axis([0 length(rf_size) 0 length(rf_size)]);
    export_fig(['./figures/NL_circle.png']); close(cf);
    
    cnt = 1;
    for s = 1:length(ssigs)                                 % for each surround scale
        
        % Gaussian rf
        surround2       = make_2D_gaussian(rf_size,csig*ssigs(s));
        srf2            = center - surround2;
        
        % plot RFs
        cin = hslcolormap(abs(round(1500*(min(srf2(:))/max(srf2(:))))),0.5,0.6,[1 0.5]); csrf2 = [flipud(cin); cex];
        sf = figure; hold on; setupfig(6,6,12); axis square; surf(srf2,'edgecolor','none');
        colormap(csrf2); set(gca,'Xtick',[],'YTick',[],'ZTick',0); axis([0 length(rf_size) 0 length(rf_size)]);
        export_fig(['./figures/filter_rat' num2str(ssigs(s)) '.png']); close(sf);
        
        % open figures for plotting CRFs
        pf = figure(20); hold on; setupfig(6,6,12); box on; title(['Surround/center = ' num2str(ssigs(s))]);
        hold off;
        
        % open figure for plotting max response ratios
        bf = figure(21); hold on; setupfig(6,6,12); box on; title(['Surround/center = ' num2str(ssigs(s))]);
        hold off;
        
        for a = 1:length(assigs)        % for each adapation scale
            
            % make adaptation filter
            ada2 = make_2D_gaussian(rf_size,csig*assigs(a));
            
            parfor k = 1:length(conts)     % for each contrast level
                
                % stimulus
                sON         = (gray.*conts(k).*c_mask) + gray;
                sOFF        = (gray.*-conts(k).*c_mask) + gray;
                
                % model response
                rONb(k)     = model_response('ON',sON,srf2,ada2,1);
                rOFFb(k)    = model_response('OFF',sOFF,srf2,ada2,1);
                
            end

            if(0) % add common nonlinearity
                rONb = rONb.^2;
                rOFFb = rOFFb.^2;
            end
            
            % NL indices
            pts = [0.2:0.1:0.9];
            
            for pt = 1:length(pts)
                
                % CRF slope around specified point (robust to floating pt errors)
                ONslope2    = rONb(abs(conts - (pts(pt)+0.05)) < 1e-10) - rONb(abs(conts - (pts(pt)-0.05)) < 1e-10);
                OFFslope2   = rOFFb(abs(conts - (pts(pt)+0.05)) < 1e-10) - rOFFb(abs(conts - (pts(pt)-0.05)) < 1e-10);
                
                % compute nonlinearity indices
                NLONb(cnt,pt)      = log10(ONslope2/(rONb(conts == 0.1) - rONb(conts == 0)));
                NLOFFb(cnt,pt)     = log10(OFFslope2/(rOFFb(conts == 0.1) - rOFFb(conts == 0)));
            end
    
            % compute max responses
            RmaxOFFb(cnt)   = max(rOFFb);
            RmaxONb(cnt)    = max(rONb);
            
            % plot CRFs
            figure(pf); hold on;
            fh(a) = plot(100*conts, rONb,'color',ColorIt(a));
            plot(100*conts, rOFFb,'--','color',ColorIt(a));
            ylim([0 2]);
            hold off;
            
            %plot Rmax OFF/Rmax ON
            figure(bf); hold on;
            plot(assigs(a), max(rOFFb)/max(rONb),'o','color',ColorIt(a),'markerfacecolor',ColorIt(a),'markersize',10);
            hold off;
            
            if s == 1
                % plot adaptation size
                af = figure; setupfig(6,6,12); axis square;
                figure(af); hold on; surf(ada2,'edgecolor','none');
                colormap(cada); set(gca,'Xtick',[],'YTick',[],'ZTick',[]); axis([0 length(rf_size) 0 length(rf_size)]);
                export_fig(['./figures/filter_ada' num2str(assigs(a)) '.png']); close(af);
            end
            
            cnt = cnt + 1;
            
        end
        
        figure(pf); hold on; set(gca,'ytick',[]); axis square; xlabel('Weber contrast'); ylabel('filter response'); %ylim([0 1]);
        legend(fh,cellstr(num2str(assigs', 'a=%f')));
        export_fig(['./figures/NL_CRF_rat' num2str(ssigs(s)) '.pdf']); close(pf)
        
        figure(bf); hold on; axis square; ylabel('OFF/ON at 100% contrast'); xlim([0.75 6.25]);
        export_fig(['./figures/NL_max_rat' num2str(ssigs(s)) '.pdf']); close(bf)

    end
    
    % model NL histograms
    for pt = 1:length(pts)
        figure; setupfig(6,6,20); hold on; box on; title(['DoG Filter']);
        hs1 = histogram(NLOFFb(:,pt),linspace(min(NLONb(:,pt)),max(NLOFFb(:,pt)),12),'facecolor',ColorIt('b'));
        hs2 = histogram(NLONb(:,pt),linspace(min(NLONb(:,pt)),max(NLOFFb(:,pt)),12),'facecolor',ColorIt('r'));
        
        xlim([min(NLONb(:,pt))-0.1 max(NLOFFb(:,pt))+0.1]); xlabel('NL index'); ylabel('Frequency'); axis square
        export_fig(['./figures/NL_histogram_' num2str(pts(pt)) '.pdf'],'-painters');
    end

end

%% L50 and LRFs %%%%%%%%%%%%%%%%%%%%%%%
if(1)
    
    shape   = {'high contrast spot','spot'};                    % names of stim types
    graykk  = 61;                                               % middle gray
    inON    = graykk:0.5:120;                                   % intensity for pixels ON
    inOFF   = graykk:-0.5:2;                                    % intensity for pixels OFF
    inONh   = linspace(2,120,numel(inON));                      % intensity for pixels ON (high contrast)
    inOFFh  = fliplr(linspace(2,120,numel(inON)));              % intensity for pixels OFF (high contrast)
    c_mask3 = double( abs(x) <= diam/2  & abs(y) <= diam/2 );   % square stimulus
    
    % initialize scatter plots
    offig = figure; hold on; setupfig(6,6,20); box on; title(['L50 OFF']);
    plot([0 1],[0 1],'k:'); xlabel('L50 light bg'); ylabel('L50 gray bg');
    
    onfig = figure; hold on; setupfig(6,6,20); box on; title(['L50 ON']);
    plot([0 1],[0 1],'k:'); xlabel('L50 dark bg'); ylabel('L50 gray bg');
    
    offigR = figure; hold on; setupfig(6,6,20); box on; title(['Rmax OFF']);
    plot([0 10],[0 10],'k:'); xlabel('Rmax light bg'); ylabel('Rmax gray bg');
    
    onfigR = figure; hold on; setupfig(6,6,20); box on; title(['Rmax ON']);
    plot([0 10],[0 10],'k:'); xlabel('Rmax dark bg'); ylabel('Rmax gray bg');
    
    titles = {'(A) OFF (bg = 120 cd)','(B) OFF (bg = 61 cd)','(C) ON (bg = 2 cd)','(D) ON (bg = 61 cd)'};
    
    % matrices to hold predicted naka-rushton values
    L50OFFbr    = []; L50OFFgr    = [];
    L50ONdk     = []; L50ONgr     = [];
    RmaxOFFbr   = []; RmaxOFFgr   = [];
    RmaxONdk    = []; RmaxONgr    = [];
    
    for s = 1:length(ssigs)     % for each surround scale
        
        % Gaussian rf
        surround2    = make_2D_gaussian(rf_size,csig*ssigs(s));
        srf2         = center - surround2;
        
        for a = 1:length(assigs)        % for each adaptation scale
            
            ada2  = make_2D_gaussian(rf_size,csig*assigs(a));
            
            %initialize lrf figure
            lrfig = figure; hold on; setupfig(8,8,20); box on;

            % naka rushton values for this RF
            L50ON   = []; L50OFF  = [];
            RmaxON  = []; RmaxOFF = [];
            
            cnt = 1;
            for p = 1:length(shape)     % for each stimulus
                
                stimim  = [];           % initialize matrix for images
                
                switch shape{p}                 % get stimulus and bg intensity values
                    
                    case 'high contrast spot'   % high contrast intensities
                        onvals  = inONh;    
                        offvals = inOFFh;
                        onbg    = 2;        
                        offbg   = 120;
                        
                    case 'spot'                 % normal contrast intensities
                        onvals  = inON;     
                        offvals = inOFF;
                        onbg    = graykk;   
                        offbg   = graykk;
                end
                
                parfor k = 1:length(onvals)             % for each intensity
                    
                    sON         = c_mask3;              % copy shape to ON stim
                    sOFF        = c_mask3;              % copy shape to OFF stim
                    
                    sON(c_mask3 == 1)   = onvals(k);    % fill ON shape
                    sON(c_mask3 == 0)   = onbg;         % fill ON bg
                    sOFF(c_mask3 == 1)  = offvals(k);   % fill OFF shape
                    sOFF(c_mask3 == 0)  = offbg;        % fill OFF bg
                    
                    dON(k)  = onvals(k) - onbg;         % calculate delta intensity with bg for plotting
                    dOFF(k) = offbg - offvals(k);
                    
                    rONx(k)  = model_response('ON',sON,srf2,ada2,1);       % compute filter response
                    rOFFx(k) = model_response('OFF',sOFF,srf2,ada2,1);
                    
                end
                
                
                % store images of stimulus
                for k = 1:length(inON)                  % for each intensity
                    
                    sON         = c_mask3;              % copy shape to ON stim
                    sOFF        = c_mask3;              % copy shape to OFF stim
                    
                    sON(c_mask3 == 1)   = onvals(k);    % fill ON shape
                    sON(c_mask3 == 0)   = onbg;         % fill ON bg
                    sOFF(c_mask3 == 1)  = offvals(k);   % fill OFF shape
                    sOFF(c_mask3 == 0)  = offbg;        % fill OFF bg
                    
                    %store image values
                    if ssigs(s) == 1.25 && assigs(a) == 1
                        if ~mod(inON(k),10); 
                            stimim = cat(2,stimim,cat(1,crop_edges(sOFF,350),crop_edges(sON,350))); 
                        end
                    end
                    
                end
                
                % normalized deltas
                dON     = dON./max(dON);
                dOFF    = dOFF./max(dOFF);
                
                % get L50 and Rmax for ON and OFF
                L50ON(cnt)      = dON(abs(rONx - max(rONx)/2) == min(abs(rONx - max(rONx)/2)));
                L50OFF(cnt)     = dON(abs(rOFFx - max(rOFFx)/2) == min(abs(rOFFx - max(rOFFx)/2)));
                RmaxON(cnt)     = max(rONx);
                RmaxOFF(cnt)    = max(rOFFx);
                
                % save stimulus
                if ssigs(s) == 1.25 && assigs(a) == 1
                    
                    figure; hold on; setupfig(12,4,20); axis image; colormap gray; set(gca,'ytick',[],'xtick',[]);
                    imagesc(stimim);
                    export_fig(['./figures/nakaLRFstim_' shape{p} '.png']);

                end
                
                figure(lrfig); hold on;
                
                subplot(2,2,cnt); hold on; box on;
                h(2) = plot(offvals, rOFFx,'--','color','k');
                set(gca,'Xtick',[2 20 40 60 80 100 120]);
                xlim([min(offvals) max(offvals)]); ylim([0 max(RmaxON)*1.1]); set(gca','Ytick',[]);
                set(gca, 'xdir','reverse'); axis square
                xlabel('spot luminance'); ylabel('filter response');
                
                subplot(2,2,cnt+2); hold on; box on;
                h(1) = plot(onvals, rONx,'color',ColorIt('r'));
                set(gca,'Xtick',[2 20 40 60 80 100 120]);
                xlim([min(onvals) max(onvals)]); ylim([0 max(RmaxON)*1.1]);
                set(gca','Ytick',[]); axis square
                xlabel('spot luminance'); ylabel('filter response');
                
                %end
                cnt = cnt + 1;
                
            end
            
            export_fig(['./figures/nakaLRF_s' num2str(ssigs(s)) '_a' num2str(assigs(a)) '.pdf']);
            close(lrfig);
            
            figure(offig); hold on;
            plot(L50OFF(1), L50OFF(2),'o','color','k','markerfacecolor','k','markersize',20);
            
            figure(onfig); hold on;
            plot(L50ON(1), L50ON(2),'o','color','k','markerfacecolor',ColorIt('r'),'markersize',20);
            
            figure(offigR); hold on;
            plot(RmaxOFF(1), RmaxOFF(2),'o','color','k','markerfacecolor','k','markersize',20);
            
            figure(onfigR); hold on;
            plot(RmaxON(1), RmaxON(2),'o','color','k','markerfacecolor',ColorIt('r'),'markersize',20);

            % store all L50s and Rmaxs
            L50OFFbr    = [L50OFFbr ; L50OFF(1)];
            L50OFFgr    = [L50OFFgr ; L50OFF(2)];
            L50ONdk     = [L50ONdk ; L50ON(1)];
            L50ONgr     = [L50ONgr ; L50ON(2)];
            
            RmaxOFFbr   = [RmaxOFFbr ; RmaxOFF(1)];
            RmaxOFFgr   = [RmaxOFFgr ; RmaxOFF(2)];
            RmaxONdk    = [RmaxONdk ; RmaxON(1)];
            RmaxONgr    = [RmaxONgr ; RmaxON(2)];
            
            
        end
        
    end
    
    display(['L50 OFF br/gr = ' num2str(mean(L50OFFbr./L50OFFgr))])
    display(['L50 ON dk/gr = ' num2str(mean(L50ONdk./L50ONgr))])
    display(['Rmax OFF br/gr = ' num2str(mean(RmaxOFFbr./RmaxOFFgr))])
    display(['Rmax ON dk/gr = ' num2str(mean(RmaxONdk./RmaxONgr))])
    
    figure(offig); hold on;
    export_fig(['./figures/nakaL50OFF.pdf']);
    
    figure(onfig); hold on;
    export_fig(['./figures/nakaL50ON.pdf']);
    
    figure(offigR); hold on;
    export_fig(['./figures/nakaRmaxOFF.pdf']);
    
    figure(onfigR); hold on;
    export_fig(['./figures/nakaRmaxON.pdf']);
    
end


%% plot Normann cone predictions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
    ad_vals = [10 1000];
    
    % plot luminance responses
    fa = figure; setupfig(6,6,12);

    semilogx(logspace(-1,5,1000),pr_response( logspace(-1,5,1000) , ad_vals(1) ), '-', 'color', ColorIt(1), 'Linewidth',2); hold on;
    semilogx(ad_vals(1),pr_response( ad_vals(1) , ad_vals(1) ), 'o', 'color', ColorIt(1), 'markerfacecolor', ColorIt(1),'markersize',10);

    semilogx(logspace(-1,5,1000),pr_response( logspace(-1,5,1000) , ad_vals(2) ), '--', 'color', ColorIt(2), 'Linewidth',2); hold on;
    semilogx(ad_vals(2),pr_response( ad_vals(2) , ad_vals(2) ), 'o', 'color', ColorIt(2), 'markerfacecolor', ColorIt(2),'markersize',10);

    set(gca,'YTick',[0 0.5 1]); xlim([0 100000]);
    xlabel('Light Intensity'); ylabel('Polarization');
    export_fig('./figures/pr_response1.pdf','-r300'); close(fa);
    
    % switch x axis to contrast
    fa = figure; setupfig(6,6,12);

    plot(100*(logspace(-1,5,1000)-ad_vals(1))/ad_vals(1),pr_response( logspace(-1,5,1000) , ad_vals(1) ), '-', 'color', ColorIt(1), 'Linewidth',2); hold on;
    plot(0,pr_response( ad_vals(1) , ad_vals(1) ), 'o', 'color', ColorIt(1), 'markerfacecolor', ColorIt(1),'markersize',10);

    plot(100*(logspace(-1,5,1000)-ad_vals(2))/ad_vals(2),pr_response( logspace(-1,5,1000) , ad_vals(2) ), '--', 'color', ColorIt(2), 'Linewidth',2); hold on;
    plot(0,pr_response( ad_vals(2) , ad_vals(2) ), 'o', 'color', ColorIt(2), 'markerfacecolor', ColorIt(2),'markersize',10);

    set(gca,'Xtick',[-100 -50 0 50 100],'YTick',[0 0.5 1]); xlim([-100 100]); ylim([0 1]);
    xlabel('Weber Contrast'); ylabel('Polarization');
    export_fig('./figures/pr_response2.pdf','-r300'); close(fa);
    
    %include ON/OFF rectificatiom
    fa = figure; setupfig(6,6,12);

    plot(abs(100*(logspace(log10(ad_vals(1)),5,1000)-ad_vals(1))/ad_vals(1)),abs(pr_response( logspace(log10(ad_vals(1)),5,1000) , ad_vals(1) )) - 0.5, '-', 'color', ColorIt('r'), 'Linewidth',2); hold on;
    plot(abs(100*(logspace(-1,log10(ad_vals(1)),1000)-ad_vals(1))/ad_vals(1)),0.5 - abs(pr_response( logspace(-1,log10(ad_vals(1)),1000) , ad_vals(1) )), '--', 'color', [0 0 0], 'Linewidth',2); hold on;

    set(gca,'Xtick',[0 50 100],'YTick',[0 0.5]); xlim([0 100]); ylim([0 0.5]);
    xlabel('Weber Contrast'); ylabel('Polarization');
    export_fig('./figures/pr_response3.pdf','-r300'); close(fa);
    
end

