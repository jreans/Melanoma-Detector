clear
close all

% S0: 2700 seconds with resize and 32 workers
% S1: 35000 seconds with resize and 32
%range =  [1:530]; % 104:105 is fast
range =  [1:120]; %6%41%[1:52];%[44]; % 104:105 is fast
load CurrentWorkingDir

cd Data
cd(folderDIR)
load GSdiag

try
    load Matrix_mIBs
    sprintf('loaded')
catch
end

cd ../..

%homeDIR = pwd;

% cd cropped
% dirdir = dir([folderDIR '*']);
% imtype = Dat_Type(1:3);
% cd ..


%% User inputs
%range_lms = 1:length(dirdir); %which pictures you want to run on
%errim = zeros(length(dirdir),1);  % Amanda, resetting this every time does not allow for [partial or cumulative runs
ensure_REDseg = 0; % uses greena dn blue layers to segment red in case where red can't be segmented
resample = 1;
FinalPixSize = 1000000;  % ***
Nvars = 100;         % number of coulmns to put into matrix for saving results
Nlayers = 3;        % number of spectral layers in image, 3 for RGB
Nims = length(Diagnoses);         % number of images in sequence to analyze
THESIZE = 200; % N Bins for both number of radii in half circle and number of pixels in 1 radial samp
clearbord = 1;      % if clearbord = 1, don't count moles that touch image border
n_ang = 18;         % analyze flip symmetry over 180 degrees
d_ang = 10;         % in 10 degree incriments
ThreshMelNetw = 0.7;  % The higher this number, the less restrictive in determining a melanocytic network
rangeMMM = 10; % flatten out DC lesion for AC melanocytic pattern by averaging over win = 2*rangeMMM

%% plotting choices
plotON_Thresh = 0;      % shows thresholding
plotON_GetStats = 0;    % shows edge fitting
plotON_Sym = 0;         % shows symmetry flipping boarkder rotation routine
plot_pix_dist = 0;          % shows pixel distribution for threshold finder
plotON_ShowCentroid =0;
plotON_Topo =0;        % Shows topographical map based on spectrum
plot_CoordXfer =0;     % shows coordinate transformation
plot_branches = 0;
SaveOutFigs = 0;

%% debugging tools
list_EdgeFits = 0;  % List the edge fit values and flags in command window
BrnchFrqAnalys = 0; % perform the more computationally expensive pig net locate
getvals = 1;
%use_whole = 1;      % Use entire image for automatic thresholding
options = optimset('Display','off', 'MaxIter',200);

%% Set up graphics
% if 1
%     fignums = [88 48];
%     Posn = [
%         10 431 1261 288
%         2573    814   825   586
%         ];
%     set(0,'Units','pixels')
%     for i_fig = 1:length(fignums)
%         set(figure(fignums(i_fig)),'position',[Posn(i_fig,:)])
%     end
% end


if 0
    fignums = [5 88 1 3 4 48 2  7  6 ];
    Posn = [
        1385.0    77.0    453.3    212.7 % 2.3333   45.6667  453.3333  212.6667
         9    21.7    1238.7    689.3%595.6667  358.3333  678.0000  352.6667
       437.6667   21.6667  810.6667  689.3333 % 2570.3    177.0    773.3    533.3
        8.3333   22.3333  416.0000  183.3333 % 3361.0    179.7    256.0    174.0
       8.3333  303.0000  415.3333  408.0000 % 3359.7    453.0    258.0    258.0
        2573.0    801.0    383.3    592.7
        2971.0    1023.7    435.3    370.0
        %3361.0    179.7    256.0    174.0
        731.0000   46.3333  279.3333  212.0000
        2971.0    810.3    434.0    116.7
        ];
    set(0,'Units','pixels')
    for i_fig = 1:length(fignums)
        set(figure(fignums(i_fig)),'position',[Posn(i_fig,:)])
    end
    
end


%% Set Wavelength Information for Data Spec.
LayerNames = [ % this is the whole spectral image cube, could be N wavelengths
    'R'
    'G'
    'B'];
Lambda = [ % these are the peak wavelengths of the spectral images in [nm]
    633
    532
    488];
Color_Order = [ 1 2 3 ]; % Pick display "R" "G" "B" wavelength data to DISPLAY data
clr = 'rgb';
%% load Marghoob Color Samples
cd('ColorSamples')
load('blackwb0')
load('bluegraywb0')
load('dbrownwb0')
load('lbrownwb0');
load('whitewb0');
load('redwb0');
cd ..

Matrix_Out1 = zeros(Nlayers,Nims,Nvars);
Matrix_Out2 = zeros(Nlayers,Nims,Nvars);
Matrix_Out3 = zeros(Nlayers,Nims,Nvars);

using_parfor = 1;

if using_parfor
    
    
    p = gcp;
    % cd cropped
    % dirdir_dat_all = dir([ folderDIR '*.tif']);
    % tic
    % addAttachedFiles(p,{dirdir_dat_all.name})
    % toc
    % cd ..
    addAttachedFiles(p,{'F08_fitERF2.m', 'F09_Fractal2.m', 'F10_GausTail.m', 'F11_TraceBranch.m', ...
        'F12_findPeaks2D.m','F13_fitellipse.m','F12_mIBs_Seg.m'})
    cd cropped
    dirdir_dat_all = dir([ folderDIR '*.tif']);
    tic
    addAttachedFiles(p,{dirdir_dat_all.name})
    toc
    
    
    % addAttachedFiles(p,{dirdir_dat_all.name})
    
    
    %     %p = gcp; % Start parallel pool
    %     % sparpool(28);
    %     p = gcp; % Start parallel pool
    %     addAttachedFiles(p,{'F08_fitERF2.m', 'F09_Fractal2.m', 'F10_GausTail.m', 'F11_TraceBranch.m', ...
    %         'F12_findPeaks2D.m','F13_fitellipse.m','F12_mIBs_Seg.m'})
    %     cd cropped
    %     dirdir_dat_all = dir([ folderDIR '*.tif']);
    %     cd ..
    %     addAttachedFiles(p,{dirdir_dat_all.name})
    %
    
    
else
    cd cropped
end


tic

parfor i = range
    
    %if Diagnoses(i,2) == 0 | Diagnoses(i,2) == 1 % if melanoma or nevus
        tic
        try
        TempBright = [];
        Out_Branches_KEEP = [];
        Out_Branches_KEEP2 = [];
        BranchStats_KEEP2 = [];
        meanClr= [];
        meanLen= [];
        size_mole= [];
        outDIST = [];
        sizeDist = [];
        Rtheta = [];
        RADlength = [];
        var_RAD = [];
        Asym = [];
        lastpix = [];
        center = [];
        BrPxLst = [];
        Out_Branches = [];
        
        %     n_left = Nims - (Nims/2-i)
        %     t_sofar = toc;
        %     T_av = i_donesofar/t_sofar;
        %     t_left = n_left * T_av;
        
        %disp(sprintf('Beginning Lesion %3.0f, which is %3.0f of %3.0f, to finish in %4.0f minutes', i, i_donesofar,Nims, t_left));
        
        %ID = num2str(Diagnoses(i,1));
        %name = strcat('Crop', ID, '.jpg');
        % reads image ID, if there was an error in autocropping and the
        % cropped image is not in the working directory, the try-catch
        % block will put a -1 in the array errim and at the end of the
        % parfor loop, the fourth column of Diagnoses will have a -1
        %     cd(targetDIR)
        
        if i < 10
            imnum = strcat(['00000' num2str(i)]);
        elseif  i >= 10 &  i < 100
            imnum = strcat(['0000' num2str(i)]);
        elseif  i >= 100 &  i < 1000
            imnum = strcat(['000' num2str(i)]);
        elseif  i >= 1000 &  i < 10000
            imnum = strcat(['00' num2str(i)]);
        elseif  i >= 10000 &  i < 100000
            imnum = strcat(['0' num2str(i)]);
        end
        
        temp_name = dirdir_Dat(i).name;
        temp_name(length(temp_name)-3:length(temp_name)) = [];
        
        %     if temp_name(length(temp_name)) =='.'
        %         temp_name(length(temp_name)) = []; % removes '.' when file is tiff
        %     end
        name = [folderDIR '_' imnum '_' temp_name '_Crop.tif'];
        
        datt = im2double(imread(name)).*255;
        
        if sum(sum(sum(datt))) == 0
            disp(sprintf('skipping #%5.0f, black image', i))
            stoppy
        end
        
        [Ny, Nx, Nlayers] = size(datt);
        % resizes image to be FinalPixSize
        if resample
            if Ny*Nx > FinalPixSize
                resampfactor = sqrt(Ny*Nx/FinalPixSize);
                datt = imresize(datt,1/resampfactor);
                [Ny, Nx, Nlayers] = size(datt);
            end
        end
        Xkeep1 = zeros(Nlayers,1);
        Xkeep2 = Xkeep1;
        Ykeep1 = Xkeep1;
        Ykeep2 = Xkeep1;
        sz = mean([Ny Nx]);
        MinMol = round((sz/10)^2);  % the min size of the mole should be a tenth of the FOV
        
        dat4MATH = datt;
        dat4SHOW = zeros(size(dat4MATH));
        dat4SHOW(:,:,1) = dat4MATH(:,:,Color_Order(1));
        dat4SHOW(:,:,2) = dat4MATH(:,:,Color_Order(2));
        dat4SHOW(:,:,3) = dat4MATH(:,:,Color_Order(3));
        
        figure(1);clf();
        subplot(3,3,1)
        imagesc(dat4SHOW(:,:,1)/256, [0 1])
        hold('on')
        %title(sprintf([LayerNames(1) '  -  ' num2str(Lambda(1))]),'fontsize',14)
        title('red','fontsize',14)
        axis('off')
        axis('image')
        colormap(gray)
        
        subplot(3,3,2)
        imagesc(dat4SHOW(:,:,2)/256, [0 1])
        hold('on')
        %title(sprintf([LayerNames(2) '  -  ' num2str(Lambda(2))]),'fontsize',14)
        title('green','fontsize',14)
        axis('off')
        axis('image')
        colormap(gray)
        
        subplot(3,3,3)
        imagesc(dat4SHOW(:,:,3), [0 256])
        hold('on')
        %title(sprintf([LayerNames(3) '  -  ' num2str(Lambda(3))]),'fontsize',14)
        title('blue','fontsize',14)
        axis('off')
        axis('image')
        colormap(gray)
        %colorbar
        
        subplot(3,3,4)
        imagesc(datt/256)
        hold('on')
        % title(sprintf(name),'fontsize',8)
        title('clock sweep','fontsize',14)
        axis('off')
        axis('image')
        drawnow
        
        subplot(3,3,5)
        imagesc(datt/256)
        hold('on')
        %title(sprintf([name]),'fontsize',8)
        title('original','fontsize',14)
        axis('off')
        axis('image')
        drawnow
        
        RADlengths = zeros (2*THESIZE+1,Nlayers);
        MoleArea = zeros(Nlayers,1);
        i_Asym_Keep = zeros(Nlayers,1);
        flag2seg = [0 0 0]; % used to flag situations where two lesion segments are found
        
        %% go through the images at different wavelengths
        for i_layer = 1:Nlayers % for instance an RGB image would have Nlayers = 3
            
            
            try
            
            slice = datt(:,:,i_layer);
            slice0 = slice;
            n4hist = 100;
            samp=slice;
            Tfff = 1;
            
            
            if i_layer == 4
                grab_Woods
            end
            
            %             if ~using_parfor
            %                 cd ..
            %             end
            %
            %
            %             F12_mIBs_Seg
            %
            %             if ~using_parfor
            %                 cd cropped
            %             end
            %
            
            
            
            
            samp= reshape(samp,Nx*Ny,1);
            mmm= mean(samp);
            [hist_n, hist_val] = hist(samp,n4hist);
            [hist_n]= smooth(hist_n,10);
            TTT = mmm;
            
            for i_iterate = 1:10  % Implement Otsu's thresholding method
                i_low = 0;
                LowSide = 0;
                i_high = 0;
                HighSide = 0;
                for iloop = 1:length(samp)
                    if samp(iloop) < TTT
                        i_low = i_low + 1;
                        LowSide = LowSide + samp(iloop);
                    end
                    if samp(iloop) > TTT
                        i_high = i_high + 1;
                        HighSide = HighSide + samp(iloop);
                    end
                end
                TTT = (HighSide/i_high + LowSide/i_low)/2;
            end
            
            if plot_pix_dist == 1
                if i_layer == 1
                    figure(5);clf();
                else
                    figure(5)
                end
                plot(hist_val,hist_n,  [clr(i_layer) '-'],'linewidth',2)
                hold('on')
                plot(TTT,0,[clr(i_layer) 'o'],'markerfacecolor',clr(i_layer))
                xlabel('pixel brightness','fontsize',16)
                ylabel('number of pixels','fontsize',16)
                title('thresholding pixel histogram','fontsize',16)
                set(gca,'fontsize',16)
                drawnow
            end
            
            
            mole=1-im2bw(slice/max(max(max(slice))),TTT*Tfff/max(max(max(slice))));
            
            %             if ensure_REDseg
            %                 if Diagnoses(i,8) % for images that were manually selected using the blue channel to threshold
            %                     mole=1-im2bw(slice_thresh/max(max(max(slice_thresh))),TTT/max(max(max(slice_thresh))));
            %                 end
%             end
            
            if plotON_Thresh
                figure(2);clf();
                subplot(3,3,1)
                imagesc(slice)
                axis('image')
                title('Original Image')
                axis('off')
                colormap(gray)
                subplot(3,3,2)
                imagesc(mole)
                axis('image')
                title('Threshold Applied')
                axis('off')
                colormap(gray)
                drawnow
            end
            
            seD = strel('diamond',1);
            mole = bwareaopen(mole,MinMol);
            
            if plotON_Thresh
                subplot(3,3,3)
                imagesc(mole)
                axis('image')
                title('bwareaopen')
                axis('off')
                colormap(gray)
                drawnow
            end
            
            mole = imfill(mole,'holes');
            
            if plotON_Thresh
                subplot(3,3,4)
                imagesc(mole)
                axis('image')
                title('imfill')
                axis('off')
                colormap(gray)
                drawnow
            end
            
            mole = imerode(mole,seD);
            
            if plotON_Thresh
                subplot(3,3,5)
                imagesc(mole)
                axis('image')
                title('imerode')
                axis('off')
                colormap(gray)
                drawnow
            end
            
            masked = mole.*slice;
            
            %             if plotON_Thresh
            %                 subplot(3,3,6)
            %                 imagesc(masked)
            %                 axis('image')
            %                 title('masked')
            %                 axis('off')
            %                 colormap(gray)
            %             end
            
            if clearbord
                mole = imclearborder(mole,4);
                masked = mole.*slice;
                if plotON_Thresh
                    subplot(3,3,6)
                    imagesc(masked)
                    axis('image')
                    title('masked')
                    axis('off')
                    colormap(gray)
                    drawnow
                end
            end
            
            mole = bwareaopen(mole,MinMol);
            
            if i_layer == 1
                Topo = mole;
            else
                try
                    Topo = Topo + mole; % try to add if red shifted mole already existed
                catch
                    Topo = mole;
                end
            end
            
            Outline = bwperim(mole,8);
            count = 0;
            PerimPts = zeros(2,0);
            for i_x = 1:Nx % make PerminPts a vector of just the lesion edge coordinates
                for i_y = 1:Ny
                    if Outline(i_y,i_x)
                        count = count + 1;
                        PerimPts(1,count) = i_x;
                        PerimPts(2,count) = i_y;
                    end
                end
            end
            
            slice_Illus = slice0;
            slice_Illus(Outline) = 255;
            
            
            
            [B, L] = bwboundaries(mole,'nohole');
            stats = regionprops(L,'all');
            
            
            
            if length(stats) > 1
                flag2seg(i_layer) = 1;
                disp(sprintf('WARNING found TWO lesion segments in channel #%2.0f of im #%4.0f, choosing largest segment!', i_layer, i))
                area_TempSeg = [];
                for i_seek = 1:length(stats)
                    area_TempSeg(i_seek) = stats(i_seek).Area;
                end
                [trash i_bigone] = max(area_TempSeg);
                for i_seek = 1:length(stats)
                    if i_seek ~= i_bigone
                        tempPixList = stats(i_seek).PixelList;
                        for i_kill = 1:length(tempPixList)
                            mole(tempPixList(i_kill,2),tempPixList(i_kill,1)) = 0;
                        end
                    end
                end
                [B, L] = bwboundaries(mole,'nohole');
                stats = regionprops(L,'all');
            end
            
            stringdat = double(reshape(slice,Nx*Ny,1));
            var = mean(stringdat)+3*std(stringdat);
            
            %infinitewhileloop
            %         if 0
            %
            %
            %         if length(stats) == 0
            %             all_done = 0;
            %             while ~all_done
            %                 Tfff = Tfff*0.9999
            %
            %                 samp= reshape(samp,Nx*Ny,1);
            %                 mmm= mean(samp);
            %                 [hist_n, hist_val] = hist(samp,n4hist);
            %                 [hist_n]= smooth(hist_n,10);
            %                 TTT = mmm;
            %
            %                 for i_iterate = 1:10  % Implement Otsu's thresholding method
            %                     i_low = 0;
            %                     LowSide = 0;
            %                     i_high = 0;
            %                     HighSide = 0;
            %                     for iloop = 1:length(samp)
            %                         if samp(iloop) < TTT
            %                             i_low = i_low + 1;
            %                             LowSide = LowSide + samp(iloop);
            %                         end
            %                         if samp(iloop) > TTT
            %                             i_high = i_high + 1;
            %                             HighSide = HighSide + samp(iloop);
            %                         end
            %                     end
            %                     TTT = (HighSide/i_high + LowSide/i_low)/2;
            %                 end
            %
            %                 if plot_pix_dist == 1
            %                     if i_layer == 1
            %                         figure(5);clf();
            %                     else
            %                         figure(5)
            %                     end
            %                     plot(hist_val,hist_n,  [clr(i_layer) '-'],'linewidth',2)
            %                     hold('on')
            %                     plot(TTT,0,[clr(i_layer) 'o'],'markerfacecolor',clr(i_layer))
            %                     xlabel('pixel brightness','fontsize',16)
            %                     ylabel('number of pixels','fontsize',16)
            %                     title('thresholding pixel histogram','fontsize',16)
            %                     set(gca,'fontsize',16)
            %                 end
            %
            %
            %                 mole=1-im2bw(slice/max(max(max(slice))),TTT*Tfff/max(max(max(slice))));
            %
            %                 %             if ensure_REDseg
            %                 %                 if Diagnoses(i,8) % for images that were manually selected using the blue channel to threshold
            %                 %                     mole=1-im2bw(slice_thresh/max(max(max(slice_thresh))),TTT/max(max(max(slice_thresh))));
            %                 %                 end
            %                 %             end
            %
            %                 if plotON_Thresh
            %                     figure(2);clf();
            %                     subplot(3,3,1)
            %                     imagesc(slice)
            %                     axis('image')
            %                     title('Original Image')
            %                     axis('off')
            %                     colormap(gray)
            %                     subplot(3,3,2)
            %                     imagesc(mole)
            %                     axis('image')
            %                     title('Threshold Applied')
            %                     axis('off')
            %                     colormap(gray)
            %                 end
            %
            %                 seD = strel('diamond',1);
            %                 mole = bwareaopen(mole,MinMol);
            %
            %                 if plotON_Thresh
            %                     subplot(3,3,3)
            %                     imagesc(mole)
            %                     axis('image')
            %                     title('bwareaopen')
            %                     axis('off')
            %                     colormap(gray)
            %                 end
            %
            %                 mole = imfill(mole,'holes');
            %                 %stop
            %                 if plotON_Thresh
            %                     subplot(3,3,4)
            %                     imagesc(mole)
            %                     axis('image')
            %                     title('imfill')
            %                     axis('off')
            %                     colormap(gray)
            %                 end
            %
            %                 mole = imerode(mole,seD);
            %
            %                 if plotON_Thresh
            %                     subplot(3,3,5)
            %                     imagesc(mole)
            %                     axis('image')
            %                     title('imerode')
            %                     axis('off')
            %                     colormap(gray)
            %                 end
            %
            %                 masked = mole.*slice;
            %
            %                 %             if plotON_Thresh
            %                 %                 subplot(3,3,6)
            %                 %                 imagesc(masked)
            %                 %                 axis('image')
            %                 %                 title('masked')
            %                 %                 axis('off')
            %                 %                 colormap(gray)
            %                 %             end
            %
            %                 if clearbord
            %                     mole = imclearborder(mole,4);
            %                     masked = mole.*slice;
            %                     if plotON_Thresh
            %                         subplot(3,3,6)
            %                         imagesc(masked)
            %                         axis('image')
            %                         title('masked')
            %                         axis('off')
            %                         colormap(gray)
            %                     end
            %                 end
            %
            %                 mole = bwareaopen(mole,MinMol);
            %
            %                 if i_layer == 1
            %                     Topo = mole;
            %                 else
            %                     try
            %                         Topo = Topo + mole; % try to add if red shifted mole already existed
            %                     catch
            %                         Topo = mole;
            %                     end
            %                 end
            %
            %                 Outline = bwperim(mole,8);
            %                 count = 0;
            %                 PerimPts = zeros(2,0);
            %                 for i_x = 1:Nx % make PerminPts a vector of just the lesion edge coordinates
            %                     for i_y = 1:Ny
            %                         if Outline(i_y,i_x)
            %                             count = count + 1;
            %                             PerimPts(1,count) = i_x;
            %                             PerimPts(2,count) = i_y;
            %                         end
            %                     end
            %                 end
            %
            %                 slice_Illus = slice0;
            %                 slice_Illus(Outline) = 255;
            %
            %
            %
            %                 [B, L] = bwboundaries(mole,'nohole');
            %                 stats = regionprops(L,'all');
            %
            %
            %
            %                 if length(stats) > 1
            %                     flag2seg(i_layer) = 1;
            %                     disp(sprintf('WARNING found TWO lesion segments in channel #%2.0f of im #%4.0f, choosing largest segment!', i_layer, i))
            %                     area_TempSeg = [];
            %                     for i_seek = 1:length(stats)
            %                         area_TempSeg(i_seek) = stats(i_seek).Area;
            %                     end
            %                     [trash i_bigone] = max(area_TempSeg);
            %                     for i_seek = 1:length(stats)
            %                         if i_seek ~= i_bigone
            %                             tempPixList = stats(i_seek).PixelList;
            %                             for i_kill = 1:length(tempPixList)
            %                                 mole(tempPixList(i_kill,2),tempPixList(i_kill,1)) = 0;
            %                             end
            %                         end
            %                     end
            %                     [B, L] = bwboundaries(mole,'nohole');
            %                     stats = regionprops(L,'all');
            %                 end
            %
            %                 stringdat = double(reshape(slice,Nx*Ny,1));
            %                 var = mean(stringdat)+3*std(stringdat);
            %
            %                 if  length(stats) == 1
            %                     all_done = 1;
            %                 end
            %             end
            %         end
            %
            %         end
            
            if ~using_parfor
                cd ..
            end
            
            fracDIM = F09_Fractal2(mole);
            
            if ~using_parfor
                cd cropped
            end
            
            
            
            MoleArea(i_layer) = stats.Area;
            
            if plotON_Thresh
                subplot(3,3,7)
                imagesc(mole)
                axis('image')
                title('bwareaopen2')
                axis('off')
                colormap(gray)
                
                subplot(3,3,8)
                imagesc(Outline)
                axis('image')
                title('Outline')
                axis('off')
                colormap(gray)
                
                subplot(3,3,9)
                imagesc(slice_Illus)
                axis('image')
                title('Marked')
                axis('off')
                colormap(gray)
                drawnow
            end
            
            
            
            %% ellipse fit of lesion border
            lbx = B{1,1}(:,2);
            lby = B{1,1}(:,1);
            trtr = stats.Centroid;
            cx = round(trtr(1));
            cy = round(trtr(2));
            
            a = cx-min(lbx);
            b = cy-min(lby);
            theta = 0;
            start(1) = a;
            start(2) = b;
            start(3) = theta;
            %         start(4) = round(cx);
            %         start(5) = round(cy);
            
            %         t=linspace(pi,2.*pi,floor(length(lbx)./2))';
            %         t2=linspace(0,pi,ceil(length(lbx)./2))';
            %         t = cat(1,t,t2);
            %         cd(homeDIR)
            %cd(homeDIR)
            
            if ~using_parfor
                cd ..
            end
            
            [resy, fval, exitflag] = fminsearch('F13_fitellipse', start, options, lbx,lby,cx,cy);
            
            if ~using_parfor
                cd cropped
            end
            
            a = resy(1);
            b = resy(2);
            theta = resy(3);
            %         cx = resy(4);
            %         cy = resy(5);
            
            if a>b
                eccentricity = sqrt(1-(b.^2./a.^2)); %biomarker
            else
                eccentricity = sqrt(1-(a.^2./b.^2));
            end
            ellerr = fval; % biomarker
            
            %close all
            figure(1)
            subplot(3,3,i_layer)
            %imagesc(slice./255);hold on
            plot(cx,cy,'w*');hold on
            plot(lbx,lby,'w-'); hold on
            t = linspace(0,2.*pi,1000);
            xell = a*cos(t)*cos(-theta) - b*sin(t)*sin(-theta) + cx;
            yell = a*cos(t)*sin(-theta) + b*sin(t)*cos(-theta) + cy;
            plot(xell,yell,'w-');
            
            %% histogram biomarkers
            n4hist = 50; %lowered to smooth histograms
            samp=slice;
            samp = samp./255;
            samp= reshape(samp,Nx*Ny,1);
            mmm= mean(samp);
            [hist_n, hist_val] = hist(samp,n4hist);
            [hist_n]= smooth(hist_n,10);
            
            [pks, locs] = findpeaks(hist_n,'Threshold',5, ...
                'MinPeakProminence',100,'MinPeakWidth',3);
            try
                difcrestn = mean(diff(pks)); %biomarker: mean diff between height of peaks
                maxdifcrestn = max(diff(pks)); %biomarker: max diff between height of peaks
            catch
            end
            
            numcrest = length(pks); %biomarker: number of crests
            crestloc = hist_val(locs);
            mcrestloc = mean(crestloc); %biomarker: mean loc of crests
            mincrestloc = min(crestloc); %biomarker
            try
                mdifcrestloc = mean(diff(crestloc)); %biomarker: mean gap between crests
            catch
            end
            absmiddifcrestloc = sum(abs(crestloc - median(hist_val))); %biomarker: sum absolute value gap between crests and center
            maxabsmiddifcrestloc = max(abs(crestloc - median(hist_val))); %biomarker
            if numcrest == 1
                difcrestn = 0;
                mdifcrestloc = 0;
                maxdifcrestn = 0;
            elseif numcrest == 0
                difcrestn = 0;
                maxdifcrestn = 0;
                crestloc = 0;
                mcrestloc = 0;
                mincrestloc = 0;
                mdifcrestloc = 0;
                absmiddifcrestloc = 0;
                maxabsmiddifcrestloc = 0;
            end
            
            %% analyze lesion segment pixels
            PixList = stats.PixelList;
            nn = length(PixList);
            sampled = zeros(nn,1);
            for ii = 1:nn
                sampled(ii) = slice0(PixList(ii,2),PixList(ii,1));
            end
            colorVAR = std(sampled)/mean(sampled);
            
            trtr = stats.Centroid;
            X = round(trtr(1));
            Y = round(trtr(2));
            
            %% Formerly get_PigFracDim
            inputIMG = (max(max(slice0)) + imcomplement(slice0)).*mole;% + (1-mole)*mean(slice0(Outline));
            
            % figure(67);clf();
            figure(88);clf();
            subplot(3,6,7)
            imagesc(slice0)
            colormap(gray)
            axis equal image
            axis('off')
            %colorbar
            %title(sprintf(name))
            title('Gray Scale')
            drawnow
            
            if ~using_parfor
                cd ..
            end
            
            inputIMG = F10_GausTail(inputIMG, PerimPts);
            
            
            if ~using_parfor
                cd cropped
            end
            
            subplot(3,6,8)
            imagesc(inputIMG)
            colormap(gray)
            axis equal image
            axis('off')
            %colorbar
            title(sprintf(name))
            title('Gaus Trail')
            drawnow
            
            % rangeMMM is a size that should be bigger than a pigmented network branch
            rangeMMM = round(sqrt(sum(sum(mole)))/5);
            % pull out small features: normalize by local mean
            I1_smooth = inputIMG./medfilt2(inputIMG,[rangeMMM rangeMMM]);
            
            for i_seek = 1:Nx
                for j_seek = 1:Ny
                    if I1_smooth(j_seek,i_seek) == Inf | mole(j_seek,i_seek) == 0
                        I1_smooth(j_seek,i_seek) = 0;
                    end
                end
            end
            
            subplot(3,6,9)
            imagesc(I1_smooth)
            colormap(gray)
            axis equal image
            axis('off')
            title('Smoothed')
            drawnow
            
            I2_adjust = histeq(I1_smooth); % tried imadjust and adapthisteq
            
            subplot(3,6,10)
            imagesc(I2_adjust)
            colormap(gray)
            axis equal image
            axis('off')
            drawnow
            
            if ~using_parfor
                cd ..
            end
            
            fracDIM_PN = F09_Fractal2(I2_adjust);                           % Metric
            
            if ~using_parfor
                cd cropped
            end
            
            title(sprintf('Frac. D=%2.3f',fracDIM_PN))
            
            %% Formerly get_pigment_network6
            %figure(88);clf();
            subplot(3,6,1)
            
            imagesc(slice0)
            title('Original Image')
            hold('on')
            axis equal image
            axis('off')
            slice_mute = slice0; % initialize for damping edges
            drawnow
            
            %win1 = round(rangeMMM/5);
            %win1 = round(sqrt(sum(sum(mole)))/5);
            
            % USER CHOICES  (NO GOOD)
            win1 = 7;  % to be changed if BrnchFrqAnalys
            win_samp = 35;
            
            if BrnchFrqAnalys
                % Find Spatial Frequency of Local Pigmented Network
                Num = [0.0163301089994429,0.00982883411072841,0.00826706738833368,0.00220069010574637,-0.00914536820159954,-0.0256985750561842,-0.0463664083595309,-0.0691255977767072,-0.0912608086038331,-0.109808055875720,-0.122171958645711,0.873499965440808,-0.122171958645711,-0.109808055875720,-0.0912608086038331,-0.0691255977767072,-0.0463664083595309,-0.0256985750561842,-0.00914536820159954,0.00220069010574637,0.00826706738833368,0.00982883411072841,0.0163301089994429];
                h = ftrans2(Num); % turns the 1D filter Num into an image
                spatial_frequency = zeros(size(2*win_samp+1,3),1);
                freqMAP = zeros(size(slice0)); % initialize
                freqMAP_out = zeros(Ny,Nx,3); % initialize
            end
            
            for ix = 1:Nx%win1+1:Nx-win1-1
                for iy = 1:Ny%win1+1:Ny-win1-1
                    if BrnchFrqAnalys
                        if iy-win_samp > 0 && iy+win_samp <= Ny & ix-win_samp > 0 & ix+win_samp <= Nx
                            I1 = double(slice0(iy-win_samp:iy+win_samp,ix-win_samp:ix+win_samp));
                            % Remove the mean (eliminate DC bias in Spectrum)
                            I1_filt = conv2(I1,h,'valid');
                            I1_unbiased = I1_filt - mean(I1_filt(:));
                            
                            % Compute Spectrum
                            [H,f1,f2] = freqz2(I1_unbiased,512,512);
                            psd = H.*conj(H);
                            
                            % Find Spectrum Peaks
                            
                            if ~using_parfor
                                cd ..
                            end
                            
                            
                            pk_location = F12_findPeaks2D(psd,.2);
                            
                            if ~using_parfor
                                cd cropped
                            end
                            
                            % Sort Peaks by Strongest
                            ind = sub2ind([512 512],pk_location(:,2),pk_location(:,1));
                            [vals,sortidx] = sort(psd(ind),'descend');
                            pk_location = pk_location(sortidx,:);
                            
                            % Calculate Spatial Frequencies
                            freq = zeros(size(pk_location,1),1);
                            ang = zeros(size(pk_location,1),1);
                            for ii = 1:size(pk_location,1)
                                center(1) = pk_location(ii,1)-257;
                                center(2) = 257-pk_location(ii,2);
                                ang(ii) = atan2(center(2),center(1))*180/pi;
                                
                                if center(1) == 0
                                    idx = round(center(2) + 257);
                                elseif center(2) == 0
                                    idx = round(center(1) + 257);
                                else
                                    idx = round(center(1)/cosd(ang(ii)) + 257);
                                end
                                
                                if idx >1 && idx < 512
                                    freq(ii) = abs(2./f1(idx));
                                else
                                    freq(ii) = 2;
                                end
                            end
                            [ufreq,ia,~] = unique(freq,'stable');
                            
                            % Use a voting scheme to pick the common frequency
                            [cnts,binc] = hist(ufreq);
                            SpFrq = binc(cnts == max(cnts));
                            
                            % sometimes there is a tie for first place where two
                            % positions have an equal high count, the following takes the mean in that
                            % case...
                            if length(SpFrq) > 1
                                SpFrq2 = mean(SpFrq);
                                SpFrq = [];
                                SpFrq = SpFrq2;
                                SpFrq2 = [];
                            end
                            
                            win1 = round(5*SpFrq);
                            freqMAP(iy,ix) = win1;
                        end
                        
                    end
                    
                    if win1 > win_samp
                        win1 = win_samp;
                    end
                    
                    stringdat = [];
                    if iy-win1 > 0 & iy+win1 < Ny & ix-win1 > 0 & ix+win1 < Nx
                        stringdat = reshape(slice0(iy-win1:iy+win1,ix-win1:ix+win1), (2*win1+1)^2,1);
                        i_count = 0;
                        for i_SD = 1:length(stringdat)
                            if stringdat(i_SD) > slice0(iy,ix)
                                i_count = i_count + 1;
                            end
                        end
                        slice_mute(iy,ix) = i_count;
                    end
                end
                
            end
            
            if BrnchFrqAnalys
                freqMAP_out(:,:,i_layer) = freqMAP;  % save frequency characteristics of local layer
            end
            
            working_im = slice_mute.*mole;
            Direct_To_Watershed = working_im;
            %working_im = 1-(working_im/max(max(working_im)));
            working_im = working_im + (1-mole)*max(max(working_im));
            
            subplot(3,6,2)
            imagesc(working_im)
            title('Ranked Image')
            axis equal image
            axis('off')
            drawnow
            %working_im = -(working_im-max(max(working_im)));
            working_im = medfilt2(working_im, [2 2]);
            
            subplot(3,6,3)
            imagesc(working_im)
            title('Smoothed ')
            axis equal image
            axis('off')
            drawnow
            
            % http://www.mathworks.com/company/newsletters/articles/the-watershed-transform-strategies-for-image-segmentation.html
            working_im_pos = imhmin(working_im,win1/2);
            working_im_neg = imhmin(imcomplement(working_im),win1/2);
            
            subplot(3,6,4)
            imagesc(working_im_pos)
            title('imhmin')
            axis equal image
            axis('off')
            drawnow
            
            %WSimP = watershed(Direct_To_Watershed);
            WSimP = watershed(working_im_pos);
            WaterShed_pos = double(WSimP).*mole;
            
            WSimN = watershed(working_im_neg);
            WaterShed_neg = double(WSimN).*mole;
            
            
            subplot(3,6,5)
            imagesc(WaterShed_pos)
            title('Watershed')
            axis equal image
            axis('off')
            colormap(gray)
            drawnow
            
            Pig_net_pos = zeros(size(WaterShed_pos)); % initialize
            Pig_net_neg = zeros(size(WaterShed_neg)); % initialize
            
            for ix = 1:Nx
                for iy = 1:Ny
                    if WaterShed_pos(iy,ix) == 0 & mole(iy,ix)
                        Pig_net_pos(iy,ix) = 1;
                    end
                    
                    if WaterShed_neg(iy,ix) == 0 & mole(iy,ix)
                        Pig_net_neg(iy,ix) = 1;
                    end
                end
            end
            
            subplot(3,6,6)
            imagesc(Pig_net_pos)
            title('Pigmented Network')
            axis equal image
            axis('off')
            colormap(gray)
            drawnow
            
            subplot(3,6,17)
            imagesc(slice0)
            hold('on')
            colormap(gray)
            axis equal image
            [iy,ix] = find(Pig_net_neg);
            plot(ix,iy,'g.','markersize',2)
            [iy,ix] = find(Pig_net_pos);
            plot(ix,iy,'r.','markersize',4)
            skel = Pig_net_pos;
            title('Neg. Net (GREEN)')
            drawnow
            
            
            % Charlie, skel is the input to the branch analyais
            
            mn=bwmorph(skel,'branchpoints');
            [row column] = find(mn);
            branchPts0    = [row column];
            endImg    = bwmorph(skel, 'endpoints');
            [row column] = find(endImg);
            endPts       = [row column];
            
            %% branchPts counts elbows, find and eliminate elbows...
            branch_Red_Flag = zeros(length(branchPts0),1);
            i_bad = 0;
            for i_check = 1:length(branchPts0)
                xpos = branchPts0(i_check,2);
                ypos = branchPts0(i_check,1);
                summ = skel(ypos-1,xpos-1) + skel(ypos-1,xpos) +skel(ypos-1,xpos+1) +skel(ypos,xpos+1) +skel(ypos+1,xpos+1) +skel(ypos+1,xpos) +skel(ypos+1,xpos-1) +skel(ypos,xpos-1);
                if summ == 2
                    branch_Red_Flag(i_check) = 1;
                    i_bad = i_bad+1;
                end
                
                if summ == 3
                    summ2 =  skel(ypos-1,xpos)  +skel(ypos,xpos+1)  +skel(ypos+1,xpos)  +skel(ypos,xpos-1);
                    if summ2 ~= 3
                        branch_Red_Flag(i_check) = 1;
                        i_bad = i_bad+1;
                    end
                end
                
                if summ == 4
                    summ2 =  skel(ypos-1,xpos)  +skel(ypos,xpos+1)  +skel(ypos+1,xpos)  +skel(ypos,xpos-1);
                    if summ2 ~= 4
                        branch_Red_Flag(i_check) = 1;
                        i_bad = i_bad+1;
                    end
                    if summ2 == 3
                        branch_Red_Flag(i_check) = 0;
                        i_bad = i_bad-1;
                    end
                end
                
            end
            
            branchPts = zeros(length(branchPts0)-i_bad,2);
            i_count = 0;
            
            for i_check = 1:length(branchPts0)
                if branch_Red_Flag(i_check) == 0
                    i_count = i_count + 1;
                    branchPts(i_count,1) = branchPts0(i_check,1);
                    branchPts(i_count,2) = branchPts0(i_check,2);
                end
            end
            
            %% reconstruct branch segments
            Out_Branches = [];
            skel2 = skel; % initialize
            i_branch = 0;
            
            for i_look = 1:length(branchPts)
                
                xpos = branchPts(i_look,2);
                ypos = branchPts(i_look,1);
                skel2(ypos,xpos) = 0;
                
                jam_node = 1;
                while jam_node == 1
                    summ2 =  skel2(ypos-1,xpos)  + skel2(ypos,xpos+1)  +skel2(ypos+1,xpos)  +skel2(ypos,xpos-1);
                    if sum(summ2) == 0
                        jam_node = 0;
                    end
                    if skel2(ypos-1,xpos)
                        Current_branch_length = 1;
                        skel2(ypos-1,xpos) = 0;
                        lastpix(1) = ypos-1;
                        lastpix(2) = xpos;
                        BrPxLst = [];
                        BrPxLst(Current_branch_length,1) = ypos-1;
                        BrPxLst(Current_branch_length,2) = xpos;
                        
                        if ~using_parfor
                            cd ..
                        end
                        
                        [skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts] = F11_TraceBranch(skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts);
                        
                        if ~using_parfor
                            cd cropped
                        end
                        
                        i_branch = i_branch+1;
                        Out_Branches(i_branch).pixlist = BrPxLst;
                    end
                    if skel2(ypos,xpos+1)
                        Current_branch_length = 1;
                        skel2(ypos,xpos+1) = 0;
                        lastpix(1) = ypos;
                        lastpix(2) = xpos+1;
                        BrPxLst = [];
                        BrPxLst(Current_branch_length,1) = ypos;
                        BrPxLst(Current_branch_length,2) = xpos+1;
                        
                        if ~using_parfor
                            cd ..
                        end
                        
                        [skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts] = F11_TraceBranch(skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts);
                        
                        if ~using_parfor
                            cd cropped
                        end
                        
                        i_branch = i_branch+1;
                        Out_Branches(i_branch).pixlist = BrPxLst;
                    end
                    if skel2(ypos+1,xpos)
                        Current_branch_length = 1;
                        skel2(ypos+1,xpos) = 0;
                        lastpix(1) = ypos+1;
                        lastpix(2) = xpos;
                        BrPxLst = [];
                        BrPxLst(Current_branch_length,1) = ypos+1;
                        BrPxLst(Current_branch_length,2) = xpos;
                        
                        if ~using_parfor
                            cd ..
                        end
                        
                        [skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts] = F11_TraceBranch(skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts);
                        
                        if ~using_parfor
                            cd cropped
                        end
                        
                        i_branch = i_branch+1;
                        Out_Branches(i_branch).pixlist = BrPxLst;
                    end
                    if skel2(ypos,xpos-1)
                        Current_branch_length = 1;
                        skel2(ypos,xpos-1) = 0;
                        lastpix(1) = ypos;
                        lastpix(2) = xpos-1;
                        BrPxLst = [];
                        BrPxLst(Current_branch_length,1) = ypos;
                        BrPxLst(Current_branch_length,2) = xpos-1;
                        
                        if ~using_parfor
                            cd ..
                        end
                        
                        [skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts] = F11_TraceBranch(skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts);
                        
                        if ~using_parfor
                            cd cropped
                        end
                        
                        i_branch = i_branch+1;
                        Out_Branches(i_branch).pixlist = BrPxLst;
                    end
                end
                
            end
            
            datt_mark_two = dat4SHOW./max(max(max(dat4SHOW)));  % load original image
            for i_x = 1:Nx
                for j_y = 1:Ny
                    if skel(j_y,i_x)
                        datt_mark_two(j_y,i_x,1) = 0;  %mark branch pixel
                        datt_mark_two(j_y,i_x,2) = 0;  %mark branch pixel
                        datt_mark_two(j_y,i_x,3) = 0;  %mark branch pixel
                    end
                end
            end
            
            figure(88)
            subplot(3,6,12)
            imagesc(datt_mark_two)
            hold('on')
            axis equal image
            axis('off')
            drawnow
            
            %% make BranchStats [ a b c d e f g ] for each pig. net. segment
            %clear BranchStats
            BranchStats = zeros(length(Out_Branches),7);
            IgSm = 2;  % IgSm = ignore small branches, less than IgSm pixels-long
            % = matrix to save, for each branch, the following:
            % [x_coord y_coord r_coord T-coord length width color]
            % a = center y position
            % b = center x position
            % c = radial distance from lesion centroid
            % d = angle from noon W/RS to lesion centroid
            % e = length of branch
            % f = mean branch brightness
            % g = standard dev of branch brightness
            
            for i_seg = 1:length(Out_Branches)%1:round(length(sts)/20):length(sts)
                
                %             clear CurrentPixList
                CurrentPixList = Out_Branches(i_seg).pixlist;
                if length(CurrentPixList) > IgSm % selects pigment branches that are >2 pixels long
                    BranchStats(i_seg,1) = mean(CurrentPixList(:,1));
                    BranchStats(i_seg,2) = mean(CurrentPixList(:,2));
                    Del_x = X - BranchStats(i_seg,1);
                    Del_y = Y - BranchStats(i_seg,2) ;
                    R = sqrt(Del_x^2 + Del_y^2);
                    theta = atan2(Del_y, Del_x)*180/pi;
                    BranchStats(i_seg,3) = R;
                    if (theta + 180) < 270   %% something wrong, wraparond mismatch
                        BranchStats(i_seg,4) = (theta + 180) + 90;
                    else
                        BranchStats(i_seg,4) = (theta + 180) + 90 - 360;
                    end
                    BranchStats(i_seg,5) = length(CurrentPixList);  % length = Perimeter/2 for a line
                    bright = zeros(BranchStats(i_seg,5),1);
                    i_zeros = 0;  % BUG
                    flag_zeros = 0; % BUG
                    for i_trg = 1:BranchStats(i_seg,5)
                        if CurrentPixList(i_trg,2) > 0 & CurrentPixList(i_trg,1) > 0% BUG
                            bright(i_trg) = slice0(CurrentPixList(i_trg,1),CurrentPixList(i_trg,2));
                        else % BUG
                            i_zeros = i_zeros + 1; % BUG
                            flag_zeros = 1; % BUG
                        end
                    end
                    
                    if flag_zeros == 0
                        BranchStats(i_seg,6) = mean(bright);
                        BranchStats(i_seg,7) = std(bright);
                    else
                        i_count = 0;
                        for i_trg = 1:length(bright)
                            if bright(i_trg) > 0
                                i_count = i_count + 1;
                                TempBright(i_count) = bright(i_trg);
                            end
                        end
                        BranchStats(i_seg,6) = mean(TempBright);
                        BranchStats(i_seg,7) = std(TempBright);
                    end
                    
                    if 1 == 0 % debug
                        %                     sts(i_seg).Centroid
                        %                     plot(ctrs(1),ctrs(2),'gd')
                        %                     plot([X ctrs(1)],[Y ctrs(2)],'r-')
                        %                     hold('on')
                        %                     text(ctrs(1)-10,ctrs(2), sprintf('%5.0f' , sts(i_seg).Perimeter))
                        %                     text(ctrs(1),ctrs(2), sprintf('%3.2f' , theta),'BackgroundColor',[1 1 1])
                        %                     %clear pixL
                        %                     pixL = sts(i_seg).PixelList;
                    end
                end
                
            end
            ii_a1 = 0;
            ddy = sum(sum(BranchStats,2)>0);
            Branch_Lengths = zeros(ddy,1);
            Branch_Bright_mean = zeros(ddy,1);
            Branch_Bright_var = zeros(ddy,1);
            BranchStats_KEEP = zeros(ddy,7);
            
            %filter out zeros
            for i_seek = 1:length(BranchStats)
                if BranchStats(i_seek,5) > 0
                    ii_a1 = ii_a1 + 1;
                    Out_Branches_KEEP(ii_a1).pixlist = Out_Branches(i_seek).pixlist;
                    BranchStats_KEEP(ii_a1,:) = BranchStats(i_seek,:);
                end
            end
            % varPOP_length = std(BranchStats_KEEP(:,5))/mean(BranchStats_KEEP(:,5));  % METRIC
            % varMEAN_bright = std(BranchStats_KEEP(:,6))/mean(BranchStats_KEEP(:,6)); % METRIC
            % varOVERALL_bright = mean(BranchStats_KEEP(:,7)./BranchStats_KEEP(:,6)); % METRIC
            [BrightHist X_BH] = hist(BranchStats_KEEP(:,6),100);
            BrightHist_SM = smooth(BrightHist,10);
            [peaksBH peakLOC ww] = findpeaks(BrightHist_SM);
            Branch_Level_mean = mean(BranchStats_KEEP(:,6));
            Branch_Level_std = std(BranchStats_KEEP(:,6));
            BranchFilt = zeros(length(BranchStats_KEEP),1);
            BranchStats_KEEP(:,1) = round(BranchStats_KEEP(:,1)); % round x and y centroids to nearest pixel
            BranchStats_KEEP(:,2) = round(BranchStats_KEEP(:,2));
            
            for i_look = 1:length(BranchFilt)
                if BrnchFrqAnalys
                    win2 = freqMAP(BranchStats_KEEP(i_look,1),BranchStats_KEEP(i_look,2));
                    win2a = Nx; % initialize to something big
                    win2b = Nx; % initialize to something big
                    win2c = Nx; % initialize to something big
                    win2d = Nx; % initialize to something big
                    flag_cutit = 0;
                    if BranchStats_KEEP(i_look,1) - win2 < 1
                        win2a = BranchStats_KEEP(i_look,1) - 1;
                        flag_cutit = 1;
                    end
                    if BranchStats_KEEP(i_look,2) - win2 < 1
                        win2b = BranchStats_KEEP(i_look,2) - 1;
                        flag_cutit = 1;
                    end
                    if BranchStats_KEEP(i_look,1) + win2 > Ny
                        win2c =  Ny - BranchStats_KEEP(i_look,1) -1;
                        flag_cutit = 1;
                    end
                    if BranchStats_KEEP(i_look,2) + win2 > Nx
                        win2d =  Nx -BranchStats_KEEP(i_look,2) -1;
                        flag_cutit = 1;
                    end
                    if flag_cutit
                        win2 = min([win2a win2b win2c win2d]);
                    end
                else
                    win2 = 2*win1;
                end
                
                % The following if statements catches an error where a branch
                % is too close to the border for the rolling window
                if BranchStats_KEEP(i_look,1)-win2 < 1 | BranchStats_KEEP(i_look,1) + win2 > Ny | BranchStats_KEEP(i_look,2)-win2 < 1 | BranchStats_KEEP(i_look,2) + win2 > Nx
                    win2 = min([BranchStats_KEEP(i_look,1)-1, BranchStats_KEEP(i_look,2)-1, Ny-BranchStats_KEEP(i_look,1), Nx-BranchStats_KEEP(i_look,2)]);
                end
                
                local_level = mean(mean( slice0(round(BranchStats_KEEP(i_look,1)+[-win2:win2]), round(BranchStats_KEEP(i_look,2)+[-win2:win2])) ));
                local_level_std = std(mean(slice0(round(BranchStats_KEEP(i_look,1)+[-win2:win2]), round(BranchStats_KEEP(i_look,2)+[-win2:win2])) ));
                if BranchStats_KEEP(i_look,6) > local_level
                    BranchFilt(i_look) = 1;
                    %figure(21)
                    if plot_branches
                        figure(88)
                        subplot(3,6,17)
                        plot(BranchStats_KEEP(i_look,2),BranchStats_KEEP(i_look,1),'kx')
                        drawnow
                    end
                end
            end
            
            %title('Brances (yellow), nodes (red). bright branches (black X)')
            axis equal image
            axis('off')
            datt_mark_three = dat4SHOW./max(max(max(dat4SHOW)));  % load original image
            
            %filter out bright branches
            ii_a1 = 0;
            for i_seek = 1:length(BranchStats_KEEP)
                if BranchFilt(i_seek) == 0
                    ii_a1 = ii_a1 + 1;
                    Out_Branches_KEEP2(ii_a1).pixlist = Out_Branches_KEEP(i_seek).pixlist;
                    BranchStats_KEEP2(ii_a1,:) = BranchStats_KEEP(i_seek,:);
                    for i_kill = 1:length(Out_Branches_KEEP(i_seek).pixlist)
                        if min(Out_Branches_KEEP(i_seek).pixlist(:,1)) > 0
                            datt_mark_three(Out_Branches_KEEP(i_seek).pixlist(i_kill,1),Out_Branches_KEEP(i_seek).pixlist(i_kill,2),1) = 0;  %mark branch pixel
                            datt_mark_three(Out_Branches_KEEP(i_seek).pixlist(i_kill,1),Out_Branches_KEEP(i_seek).pixlist(i_kill,2),2) = 0;  %mark branch pixel
                            datt_mark_three(Out_Branches_KEEP(i_seek).pixlist(i_kill,1),Out_Branches_KEEP(i_seek).pixlist(i_kill,2),3) = 0;  %mark branch pixel
                        end
                    end
                end
            end
            varPOP_length = std(BranchStats_KEEP2(:,5))/mean(BranchStats_KEEP2(:,5));  % METRIC
            varMEAN_bright = std(BranchStats_KEEP2(:,6))/mean(BranchStats_KEEP2(:,6)); % METRIC
            tempstr = BranchStats_KEEP2(:,7)./BranchStats_KEEP2(:,6);
            totnumNaN = sum(isnan(tempstr));
            tempstr2 = zeros(1,length(tempstr)-totnumNaN);
            c_seek = 0;
            for i_seek = 1:length(tempstr)
                if ~isnan(tempstr(i_seek))
                    c_seek = c_seek + 1;
                    tempstr2(c_seek) = tempstr(i_seek);
                end
            end
            
            varOVERALL_bright = mean(tempstr2); % METRIC
            
            i_count_BP = 0;
            i_count_EP = 0;
            for i_seek = 1:length(Out_Branches_KEEP)
                %clear temppixlist
                temppixlist = Out_Branches_KEEP(i_seek).pixlist;
                for j_seek = 1:length(temppixlist)
                    for k_seek = 1:length(branchPts)
                        if branchPts(k_seek,:) == temppixlist(j_seek,:)
                            i_count_BP = i_count_BP + 1;
                        end
                    end
                end
                for j_seek = 1:length(temppixlist)
                    for k_seek = 1:length(endPts)
                        if endPts(k_seek,:) == temppixlist(j_seek,:)
                            i_count_EP = i_count_EP + 1;
                        end
                    end
                end
            end
            
            Connectedness = i_count_BP/i_count_EP;  % METRIC
            
            %figure(44)
            subplot(3,6,11)
            imagesc(datt_mark_three)
            axis equal image
            axis('off')
            drawnow
            
            [wwa wwb] = max(ww);
            i_max_peak = wwb;
            ww(wwb) = 0;
            [wwa wwb] = max(ww);
            i_2nd_peak = wwb;
            
            %figure(49) ;clf();
            figure(48) ;clf();
            subplot(6,1,6)
            plot(X_BH,BrightHist,'k-')
            hold('on')
            plot(X_BH,BrightHist_SM,'b-')
            xlabel('Branch Color')
            ylabel('#branches')
            plot(X_BH(peakLOC), peaksBH,'ro','markerfacecolor','r')
            TTTbh = mean(X_BH);
            drawnow
            
            for i_iterate = 1:10  % Implement Otsu's thresholding method
                i_low = 0;
                LowSide = 0;
                i_high = 0;
                HighSide = 0;
                for iloop = 1:length(BrightHist)
                    if BrightHist(iloop) < TTTbh
                        i_low = i_low + 1;
                        LowSide = LowSide + BrightHist(iloop);
                    end
                    if BrightHist(iloop) > TTTbh
                        i_high = i_high + 1;
                        HighSide = HighSide + BrightHist(iloop);
                    end
                end
                TTTbh = (HighSide/i_high + LowSide/i_low)/2;
            end
            
            clrs = zeros(360,1);
            lens = zeros(360,1);
            [LBS trash] = size(BranchStats_KEEP2);
            
            for i_seek = 1:LBS;
                if round(BranchStats_KEEP2(i_seek,4)+1) < 361
                    clrs(round(BranchStats_KEEP2(i_seek,4)+1)) = clrs(round(BranchStats_KEEP2(i_seek,4)+1)) + 1;
                    lens(round(BranchStats_KEEP2(i_seek,4)+1)) = lens(round(BranchStats_KEEP2(i_seek,4)+1)) + 1;
                end
            end
            
            MaxClrs = max(clrs);
            MaxLens = max(lens);
            ClrVctr = zeros(360,MaxClrs);
            LenVctr = zeros(360,MaxLens);
            clrs = zeros(360,1); % REPEATED FROM ABOVE!!!!
            lens = zeros(360,1);
            
            for i_seek = 1:LBS;
                if round(BranchStats_KEEP2(i_seek,4)+1) < 361
                    clrs(round(BranchStats_KEEP2(i_seek,4)+1)) = clrs(round(BranchStats_KEEP2(i_seek,4)+1)) + 1;
                    lens(round(BranchStats_KEEP2(i_seek,4)+1)) = lens(round(BranchStats_KEEP2(i_seek,4)+1)) + 1;
                    ClrVctr( round(BranchStats_KEEP2(i_seek,4)+1), clrs(round(BranchStats_KEEP2(i_seek,4)+1)) ) = BranchStats_KEEP2(i_seek,7);
                    LenVctr( round(BranchStats_KEEP2(i_seek,4)+1), lens(round(BranchStats_KEEP2(i_seek,4)+1)) )= BranchStats_KEEP2(i_seek,5);
                end
            end
            StdClr = zeros(1,360);
            StdLen = zeros(1,360);
            for i_seek = 1:360
                %clear inds
                inds = find(ClrVctr(i_seek,:));
                meanClr(i_seek) = mean(ClrVctr(i_seek,inds));
                StdClr(i_seek) = std(ClrVctr(i_seek,inds));
                meanLen(i_seek) = mean(LenVctr(i_seek,inds));
                StdLen(i_seek) = std(LenVctr(i_seek,inds));
            end
            
            figure(6);
            subplot(1,3,1)
            imagesc(dat4SHOW./max(max(max(dat4SHOW))))
            axis equal image
            axis('off')
            
            figure(88)
            subplot(3,6,18)
            imagesc(-skel/2, [-1 0]);
            hold('on')
            colormap(gray)
            axis equal image
            axis('off')
            plot(branchPts(:,2),branchPts(:,1),'r.');
            plot(endPts(:,2),endPts(:,1),'b.');
            axis equal
            numBRANCHES = zeros(365);
            clrBRANCHES = zeros(365);
            drawnow
            
            MaxNetY = max(endPts(:,1));
            MinNetY = min(endPts(:,1));
            MaxNetX = max(endPts(:,2));
            MinNetX = min(endPts(:,2));
            
            for i_seg = 1:length(Out_Branches_KEEP2)
                ctrXXX = mean(Out_Branches_KEEP2(i_seg).pixlist(:,2));
                text(mean(Out_Branches_KEEP2(i_seg).pixlist(:,2)),mean(Out_Branches_KEEP2(i_seg).pixlist(:,1)), sprintf('%5.0f' , length(Out_Branches_KEEP2(i_seg).pixlist)),'fontsize',4)
            end
            
            XaxVect = [1:360]; % degrees
            
            axis([MinNetX-1 MaxNetX+1 MinNetY-1 MaxNetY+1])
            
            
            figure(48);
            subplot(6,1,4)
            plot(BranchStats_KEEP2(:,4),BranchStats_KEEP2(:,5),'k*')
            hold('on')
            plot(XaxVect, smooth(meanLen,10), 'bo-','markersize',3,'markerfacecolor','b')
            % for i_seek = 1:360
            %     plot([XaxVect(i_seek) XaxVect(i_seek)], [meanLen(i_seek)-StdLen(i_seek) meanLen(i_seek)+StdLen(i_seek) ], 'b-','linewidth',3)
            % end
            % title('branch lengths')
            %xlabel('angle [degrees]')
            ylabel('length')
            drawnow
            
            subplot(6,1,5)
            plot(BranchStats_KEEP2(:,4),BranchStats_KEEP2(:,7),'k*')
            hold('on')
            plot(XaxVect, smooth(meanClr,10), 'bo-','markersize',3,'markerfacecolor','b')
            % for i_seek = 1:360
            %     plot([XaxVect(i_seek) XaxVect(i_seek)], [meanClr(i_seek)-StdClr(i_seek) meanClr(i_seek)+StdClr(i_seek) ], 'b-','linewidth',3)
            % end
            %title('branch brightnes')
            %xlabel('angle [degrees]')
            ylabel('brightnes')
            drawnow
            %axis off
            
            subplot(6,1,3)
            plot(1:length(lens),lens,'k-','linewidth',5)
            % for i_seek = 1:360
            %     plot([XaxVect(i_seek) XaxVect(i_seek)], [meanClr(i_seek)-StdClr(i_seek) meanClr(i_seek)+StdClr(i_seek) ], 'b-','linewidth',3)
            % end
            %title('Number of Branches')
            ylabel('#branches')
            %ylabel('brightnes')
            %axis off
            drawnow
            
            rangeClr = (max(smooth(meanClr,10)) - min(smooth(meanClr,10))) / mean(smooth(meanClr,10));
            rangeLen = (max(smooth(meanLen,10)) - min(smooth(meanLen,10))) / mean(smooth(meanLen,10));
            rangeNum = (max(smooth(lens,10)) - min(smooth(lens,10))) / mean(smooth(lens,10));
            
            
            %         %         %% set up figures for showing pigmented network
            %                 fignums = [88 46 45 44 49 67 21 47 48];
            %                 Posn = [64   669   560   420
            %                     648   673   560   420
            %                     648   673   560   420
            %                     648   673   560   420
            %                     1240         674         560         420
            %                     85    89   560   420
            %                     672    88   560   420
            %                     2492          19         985         790
            %                     1263          86         533         421
            %                     ];
            %                 set(0,'Units','pixels')
            %                 for i_fig = 1:length(fignums)
            %                     set(figure(fignums(i_fig)),'position',[Posn(i_fig,:)])
            %                 end
            %         %
            Just_Mole = masked.*mole;
            size_mole(i_layer) = sum(sum(mole));
            
            if plotON_ShowCentroid
                %figure(8);clf();
                figure(88)
                subplot(3,6,13)
                imagesc(Just_Mole)
                axis equal
                axis('off')
                colormap(gray)
                %colorbar
                title('original')
                drawnow
            end
            
            BWw = Just_Mole > 0 ;
            minJM = min(min(Just_Mole));
            Just_Mole = Just_Mole - minJM;
            Just_Mole = Just_Mole.*mole;
            
            if plotON_ShowCentroid
                figure(88)
                subplot(3,6,14)
                imagesc(Just_Mole)
                axis equal
                axis('off')
                colormap(gray)
                %colorbar
                title('zeroed out')
                drawnow
            end
            
            Just_Mole = Just_Mole/max(max(Just_Mole)); % Normalize
            
            if plotON_ShowCentroid
                figure(88)
                subplot(3,6,15)
                imagesc(Just_Mole)
                axis equal 
                axis('off')
                colormap(gray)
                %colorbar
                title('Normalized')
                drawnow
            end
            
            Just_Mole = 1-Just_Mole; % Invert
            Just_Mole = Just_Mole.*mole;
            if plotON_ShowCentroid
                figure(88)
                subplot(3,6,16)
                imagesc(Just_Mole)
                hold('on')
                axis equal image
                axis('off')
                colormap(gray)
                title('Inverted')
                drawnow
            end
            
            statsWeighted = regionprops(BWw, Just_Mole, {'Centroid','WeightedCentroid'});
            tempCTR  =   statsWeighted.WeightedCentroid;
            
            Xw = round(tempCTR(1));
            Yw = round(tempCTR(2));
            if i_layer == 1
                centers_spectral = zeros(Nlayers,2);
            end
            centers_spectral(i_layer,1) = Yw;
            centers_spectral(i_layer,2) = Xw;
            
            
            %         if plotON_ShowCentroid
            %             figure(8)
            %             subplot(2,2,4)
            %             plot(Xw,Yw,[clr(i_layer) '*'])
            %             drawnow
            %         end
            
            Xkeep1(i_layer) = X;
            Ykeep1(i_layer) = Y;
            Xkeep2(i_layer) = Xw;
            Ykeep2(i_layer) = Yw;
            
            sizelist2 = sort([Xw Yw]);
            nnn = sizelist2(1)-1;
            brd = stats(1).Perimeter/sqrt(stats(1).Area) - 2*pi/sqrt(pi);  % output
            
            %clear dif dif2 Assym2
            
            XXX = zeros(n_ang,1); % initialize arrays
            YYY = XXX;
            Assym2 = XXX;
            
            for ii = 1:n_ang
                deg_rot = ii*d_ang;
                B = [];
                L = [];
                rotated = [];
                ctr = [];
                %stats.Centroid = [];
                flipout2 = [];
                dif2 = [];
                rotated = logical(imrotate(mole,deg_rot,'nearest','loose')); % Nearest neighbor interpolation
                [Ny, Nx] = size(rotated);
                rotated = bwareaopen(rotated,MinMol);
                [B, L] = bwboundaries(rotated,'nohole');
                stats2 = regionprops(L,'all');
                trtr = stats2.Centroid;
                XX = round(trtr(1));
                YY = round(trtr(2));
                XXX(ii) = XX;
                YYY(ii) = YY;
                flipout2 = rotated';
                [BB, LL] = bwboundaries(flipout2,'nohole');
                stats3 = regionprops(L,'all');
                XXf = round(stats3.Centroid(1));
                YYf = round(stats3.Centroid(2));
                sizelist2 = sort([XX YY]);
                nnn = sizelist2(1)-1;
                factorBIG = 4;
                dif2 = zeros(factorBIG*nnn,factorBIG*nnn);
                
                for iii = 1:factorBIG*nnn
                    for j = 1:factorBIG*nnn
                        if YY-XXf+iii > 0 && XX-YYf+j > 0 && XX-YYf+j < Nx && YY-XXf+iii < Ny && j < Ny && iii < Nx
                            dif2(j,iii) = abs(rotated(j,iii) - flipout2(XX-YYf+j,YY-XXf+iii));
                        end
                    end
                end
                [NdiffY, NdiffX] = size(dif2);
                Assym2(ii) = sum(reshape(dif2,NdiffX*NdiffY,1))/nn;
                
                if plotON_Sym == 1
                    
                    YdiffSamp = sum(dif2,2);
                    XdiffSamp = sum(dif2,1);
                    
                    
                    %                 for iii = round(length(YdiffSamp)/2):length(YdiffSamp)
                    %                     if YdiffSamp(iii)>0
                    %                         i_Yhigh = iii;
                    %                     end
                    %                     if XdiffSamp(iii)>0
                    %                         i_Xhigh = iii;
                    %                     end
                    %                     if YdiffSamp(iii-round(length(YdiffSamp)/2)+1)>0
                    %                         i_Ylow = iii-round(length(YdiffSamp)/2)+1;
                    %                     end
                    %                     if XdiffSamp(iii-round(length(YdiffSamp)/2)+1)>0
                    %                         i_Xlow = iii-round(length(YdiffSamp)/2)+1;
                    %                     end
                    %                 end
                    
                    iii = length(YdiffSamp);
                    while YdiffSamp(iii) == 0
                        iii = iii-1;
                    end
                    i_Yhigh = iii;
                    
                    iii = 1;
                    while YdiffSamp(iii) == 0
                        iii = iii+1;
                    end
                    i_Ylow = iii;
                    
                    iii = length(YdiffSamp);
                    while XdiffSamp(iii) == 0
                        iii = iii-1;
                    end
                    i_Xhigh = iii;
                    
                    iii = 1;
                    while XdiffSamp(iii) == 0
                        iii = iii+1;
                    end
                    i_Xlow = iii;
                    
                    figure(3)
                    if ii == 1
                        figure(3)
                        clf();
                    end
                    subplot(3,6,ii)
                    imagesc(dif2)
                    hold('on')
                    axis([min(i_Xlow, i_Ylow) max(i_Xhigh,  i_Yhigh) min(i_Xlow, i_Ylow) max(i_Xhigh,  i_Yhigh)])
                    axis('off')                   
                    colormap(gray)
                    title(sprintf('A = %0.4f',Assym2(ii)),'fontsize',4)
                    plot(XX,YY,'gx')
                    drawnow
                end
            end
            
            [Big_Asym, garbage] = max(Assym2);
            [sym, i_sym] = min(Assym2);
            if i_sym == 9
                i_sym = 8;
            end
            if i_sym == 18
                i_sym = 17;
            end
            if plotON_Sym == 1
                subplot(3,6,i_sym)
                plot(XXX(i_sym),YYY(i_sym),'bo','markerfacecolor','b')
            end
            n_shift = round(90/d_ang);
            i_Asym = i_sym + n_shift;
            if i_sym > n_ang/2
                i_Asym = i_sym - n_shift;
            end
            Asym(i) = Assym2(i_Asym);
            if plotON_Sym == 1
                subplot(3,6,i_Asym)
                plot(XXX(i_Asym),YYY(i_Asym),'ro','markerfacecolor','r')
            end
            i_Asym_Keep(i_layer) = i_Asym;
            [Nyy, Nxx] = size(slice);
            ThetaTS = (i_sym*d_ang)*pi/180;
            ThetaTS_asym = (i_Asym*d_ang)*pi/180;
            
            drawnow
            
            for ix = 1:X
                xplot = X+ix;
                xplotN = X-ix;
                yp = Y-ix*tan(ThetaTS);
                yn = Y+ix*tan(ThetaTS);
                yyp = Y-ix*tan(ThetaTS_asym);
                yyn = Y+ix*tan(ThetaTS_asym);
                if round(xplot) > 0 && round(xplot) < Nxx && round(yp) > 0 && round(yp) < Nyy
                    if  mole(round(yp),round(xplot))
                        x1 = xplot;
                        y1 = yp;
                    end
                end
                if round(xplotN) > 0 && round(xplotN) < Nxx && round(yn) > 0 && round(yn) < Nyy
                    if  mole(round(yn),round(xplotN))
                        x2 = xplotN;
                        y2 = yn;
                    end
                end
                if round(xplot) > 0 && round(xplot) < Nxx && round(yyp) > 0 && round(yyp) < Nyy
                    if  mole(round(yyp),round(xplot))
                        x1_asym = xplot;
                        y1_asym = yyp;
                    end
                end
                if round(xplotN) > 0 && round(xplotN) < Nxx && round(yyn) > 0 && round(yyn) < Nyy
                    if  mole(round(yyn),round(xplotN))
                        x2_asym = xplotN;
                        y2_asym = yyn;
                    end
                end
            end
            
            diampix1 = sqrt((x1_asym-x2_asym)^2+(y1_asym-y2_asym)^2);
            diampix2 = sqrt((x1-x2)^2+(y1-y2)^2);
            diampix = (diampix1 + diampix2)/2;
            
            if i_layer == length(LayerNames)
                i_count_pixdist = 0;
                for i_layer_look = 1:length(LayerNames)
                    for j_layer_look = 1:length(LayerNames)
                        if abs(i_layer_look - j_layer_look) > 0
                            i_count_pixdist = i_count_pixdist + 1;
                            Xdist = centers_spectral(i_layer_look,2) - centers_spectral(j_layer_look,2);
                            Ydist = centers_spectral(i_layer_look,1) - centers_spectral(j_layer_look,1);
                            outDIST(i_count_pixdist) = sqrt( Xdist^2 + Ydist^2);
                            sizeDist(i_count_pixdist) = abs(size_mole(i_layer_look) - size_mole(j_layer_look));
                        end
                    end
                end
                clroff = mean(outDIST)/sqrt(mean(size_mole)); % metric based on geometric centroids
                sz_dif = mean(sizeDist)/mean(size_mole); % metric!
                sz_dif_max = max(sizeDist/mean(size_mole)); % metric!\
                %stop
            end
            
            %% formerly getSTATS  ... initializations:
            i_CountRadials = 0; % initialize for forward and backward for loops below
            XXXX = zeros(THESIZE,1);
            YYYY = XXXX;
            XXXXX = XXXX;
            YYYYY = XXXX;
            sampy = XXXX;
            dummy = XXXX;
            RADmean = zeros(THESIZE,1);
            RADstd = zeros(THESIZE,1);
            Out_slope = zeros(THESIZE,1);
            Out_slopelinmid = zeros(THESIZE,1);
            Out_slopelinbeg = zeros(THESIZE,1);
            Out_slopelinend = zeros(THESIZE,1);
            
            if plot_CoordXfer
                OutXfer = zeros(THESIZE,THESIZE); % this is the matrix for angular to cartesian transfer output
            end
            
            if getvals
                Out_flagVAL = zeros(THESIZE,2);
                Out_flagVALlinbeg = zeros(THESIZE,2);
                Out_flagVALlinmid = zeros(THESIZE,2);
                Out_flagVALlinend = zeros(THESIZE,2);
            end
            
            for theta = -pi+pi/THESIZE:pi/THESIZE:pi+pi/THESIZE
                i_CountRadials = i_CountRadials+1;
                
                for ix = 1:2*X
                    xplot = X - sin(theta)*ix;
                    yp = Y + cos(theta)*ix;
                    if round(xplot) > 0 && round(xplot) < Nxx && round(yp) > 0 && round(yp) < Nyy
                        if  mole(round(yp),round(xplot))
                            x1 = xplot;
                            y1 = yp;
                        end
                    end
                end
                
                Rtheta(i_CountRadials) = sqrt( (x1-X)^2 + (y1-Y)^2 );
                
                if plotON_GetStats == 1
                    figure(1)
                    subplot(3,3,4)
                    plot(x1,y1, [clr(i_layer) '.'],'markersize',3)
                    if i_CountRadials == 1
                        plot(X,Y,[clr(i_layer) '*'])
                        plot(Xw,Yw,[clr(i_layer) 'o'])
                    end
                    drawnow
                    pause(0.001)
                end
                
                delX = x1-X;
                delY = y1-Y;
                XXXX = round((X:((delX)/THESIZE):x1)); % for pixels in lesion
                YYYY = round((Y:((delY)/THESIZE):y1));
                XXXXX = round((X+delX/2:(delX/THESIZE):X+delX*3/2)); % for edge
                YYYYY = round((Y+delY/2:(delY/THESIZE):Y+delY*3/2));
                
                if abs(delX) < 0.1 % if the radial is straight in x direction
                    XXXX = zeros(length(YYYY),1);
                    XXXX = XXXX + X;
                    XXXXX = zeros(length(YYYYY),1);
                    XXXXX = XXXXX + X;
                end
                
                if abs(delY) < 0.1 % if the radial is straight in x direction
                    YYYY = zeros(length(XXXX),1);
                    YYYY = YYYY + Y;
                    YYYYY = zeros(length(XXXXX),1);
                    YYYYY = YYYYY + Y;
                end
                
                rngY = max(YYYYY)-min(YYYYY);
                rngX = max(XXXXX)-min(XXXXX);
                norm_size = sqrt(rngY^2+rngX^2);
                
                for i_samp = 1:THESIZE
                    sampy(i_samp) = slice(YYYY(i_samp),XXXX(i_samp));
                    if YYYYY(i_samp) > 0 && XXXXX(i_samp) > 0 && YYYYY(i_samp) < Ny && XXXXX(i_samp) < Nx
                        dummy(i_samp) = slice(YYYYY(i_samp),XXXXX(i_samp));
                    end
                end
                
                mid_dummy = min(dummy) + (max(dummy)-min(dummy))/2;
                i_middle = 0;
                
                for i_dmy = 1:length(dummy)
                    if dummy(i_dmy) < mid_dummy % find 1/2max: initial fitting param
                        i_middle = i_dmy;
                    end
                    
                    if dummy(i_dmy) < mid_dummy*1.5 % find 3/4max: initial fitting param
                        i_high = i_dmy;
                    end
                end
                
                if max(dummy) > 0
                    delta_r = dummy(i_high) - dummy(i_middle);
                    bbb = delta_r;
                    offRr = min(dummy);
                    Cmax = max(dummy);
                    bbb = -.03;
                    offRr = 100;
                    Cmax = .01;
                    
                    offRrlin = min(dummy);
                    Cmaxlin = delta_r;
                    
                    Cmaxlinbeg = dummy(length(dummy)./5) - dummy(1);
                    
                    offRrlinend = dummy(length(dummy).*4./5);
                    Cmaxlinend = dummy(end) - dummy(length(dummy).*4./5);
                    
                    
                    if dummy(round(length(dummy)/2)) > 0
                        start = [bbb offRr Cmax];
                        
                        if ~using_parfor
                            cd ..
                        end
                        
                        [resy, fval, exitflag, outMSG] = fminsearch('F08_fitERF2', start, options, dummy, THESIZE);
                        
                        if ~using_parfor
                            cd cropped
                        end
                        
                        Out_flagVAL(i_CountRadials,1) = exitflag;
                        Out_flagVAL(i_CountRadials,2) = fval;
                        bbb = resy(1);
                        offRr = resy(2);
                        Cmax = resy(3);
                        Out_slope(i_CountRadials) = bbb;   % metric transformed
                        
                        %                         start = [offRrlin Cmaxlin];
                        %                         [resy, fval, exitflag, outMSG] = fminsearch('fitERFlinmid', start, options, dummy, THESIZE);
                        %                         Out_flagVALlinmid(i_CountRadials,1) = exitflag;
                        %                         Out_flagVALlinmid(i_CountRadials,2) = fval;
                        %                         offRrlinmid = resy(1);
                        %                         Cmaxlinmid = resy(2);
                        %                         Out_slopelinmid(i_CountRadials) = Cmaxlinmid;   % metric transformed
                        
                        %                         start = [offRrlin Cmaxlinbeg];
                        %                         [resy, fval, exitflag, outMSG] = fminsearch('fitERFlinbeg', start, options, dummy, THESIZE);
                        %                         Out_flagVALlinbeg(i_CountRadials,1) = exitflag;
                        %                         Out_flagVALlinbeg(i_CountRadials,2) = fval;
                        %                         offRrlinbeg = resy(1);
                        %                         Cmaxlinbeg = resy(2);
                        %                         Out_slopelinbeg(i_CountRadials) = Cmaxlinbeg;   % metric transformed
                        %
                        %                         start = [offRrlinend Cmaxlinend];
                        %                         [resy, fval, exitflag, outMSG] = fminsearch('fitERFlinend', start, options, dummy, THESIZE);
                        %                         Out_flagVALlinend(i_CountRadials,1) = exitflag;
                        %                         Out_flagVALlinend(i_CountRadials,2) = fval;
                        %                         offRrlinend = resy(1);
                        %                         Cmaxlinend = resy(2);
                        %                         Out_slopelinend(i_CountRadials) = Cmaxlinend;   % metric transformed
                        
                        
                        if plotON_GetStats == 1
                            
                            if PCscreen_remote
                                set(figure(4),'position',1000*[ 0.0103    0.0363    1.2607    0.6747])
                            end
                            
                            figure(4)
                            subplot(2,1,2);
                            hold off
                            plot(dummy,'kx')
                            hold('on')
                            xxtemp = (1:length(dummy));
                            %ppp_erf1 = erf(( xxtemp  -round(THESIZE/2) )/bbb);
                            ppp_erf = offRr + 1./(Cmax+exp(bbb.*xxtemp));
                            plot(ppp_erf,'k-','linewidth',2)
                            title(sprintf('Lesion edge slope = %3.3f',bbb));
                            ylabel('Brightness');
                            xlabel('Pixels across lesion EDGE','fontsize',16);
                            axis([0 THESIZE min(dummy) max(dummy)])
                            drawnow
                            pause(0.001)
                        end
                    end
                    
                end
                
                if plotON_GetStats == 1
                    figure(4)
                    subplot(2,1,1)
                    hold off
                    plot(sampy,'k-','linewidth',2)
                    hold('on')
                    axis([0 THESIZE 0 256])
                    title('inside lesion')
                    ylabel('Brightness');
                    drawnow
                    pause(0.001)
                end
                RADmean(i_CountRadials) = mean(sampy);
                RADstd(i_CountRadials) = std(sampy);
                RADlength(i_CountRadials) = sqrt((x1-X)^2 + (y1-Y)^2);
                
                if plot_CoordXfer
                    OutXfer(i_CountRadials,:) = sampy;
                end
                
            end
            
            n_put = sum(Out_flagVAL(:,1));
            Sgv = zeros(n_put,1);
            Sgf = zeros(n_put,1);
            Sbv = zeros(length(Out_flagVAL)-n_put,1);
            Sbf = zeros(length(Out_flagVAL)-n_put,1);
            Out_slope2 = zeros(n_put,1);
            
            i_good = 0;
            i_bad = 0;
            for tttt = 1:length(Out_flagVAL)
                if Out_flagVAL(tttt,1)
                    i_good = i_good + 1;
                    Sgv(i_good) = Out_slope(tttt);
                    Sgf(i_good) = Out_flagVAL(tttt,2);
                end
                
                if Out_flagVAL(tttt,1) == 0
                    i_bad = i_bad + 1;
                    Sbv(i_bad) = Out_slope(tttt);
                    Sbf(i_bad) = Out_flagVAL(tttt,2);
                end
                
            end
            
            %             n_put = sum(Out_flagVALlinmid(:,1));
            %             Sgvlinmid = zeros(n_put,1);
            %             Sgflinmid = zeros(n_put,1);
            %             Sbvlinmid = zeros(length(Out_flagVALlinmid)-n_put,1);
            %             Sbflinmid = zeros(length(Out_flagVALlinmid)-n_put,1);
            %             Out_slope2linmid = zeros(n_put,1);
            %
            %             i_good = 0;
            %             i_bad = 0;
            %             for tttt = 1:length(Out_flagVALlinmid)
            %                 if Out_flagVALlinmid(tttt,1)
            %                     i_good = i_good + 1;
            %                     Sgvlinmid(i_good) = Out_slopelinmid(tttt);
            %                     Sgflinmid(i_good) = Out_flagVALlinmid(tttt,2);
            %                 end
            %
            %             end
            %
            %             n_put = sum(Out_flagVALlinbeg(:,1));
            %             Sgvlinbeg = zeros(n_put,1);
            %             Sgflinbeg = zeros(n_put,1);
            %             Sbvlinbeg = zeros(length(Out_flagVALlinbeg)-n_put,1);
            %             Sbflinbeg = zeros(length(Out_flagVALlinbeg)-n_put,1);
            %             Out_slope2linbeg = zeros(n_put,1);
            %
            %             i_good = 0;
            %             i_bad = 0;
            %             for tttt = 1:length(Out_flagVALlinbeg)
            %                 if Out_flagVALlinbeg(tttt,1)
            %                     i_good = i_good + 1;
            %                     Sgvlinbeg(i_good) = Out_slopelinbeg(tttt);
            %                     Sgflinbeg(i_good) = Out_flagVALlinbeg(tttt,2);
            %                 end
            %
            %             end
            %
            %             n_put = sum(Out_flagVALlinend(:,1));
            %             Sgvlinend = zeros(n_put,1);
            %             Sgflinend = zeros(n_put,1);
            %             Sbvlinend = zeros(length(Out_flagVALlinend)-n_put,1);
            %             Sbflinend = zeros(length(Out_flagVALlinend)-n_put,1);
            %             Out_slope2linend = zeros(n_put,1);
            %
            %             i_good = 0;
            %             i_bad = 0;
            %             for tttt = 1:length(Out_flagVALlinend)
            %                 if Out_flagVALlinend(tttt,1)
            %                     i_good = i_good + 1;
            %                     Sgvlinend(i_good) = Out_slopelinend(tttt);
            %                     Sgflinend(i_good) = Out_flagVALlinend(tttt,2);
            %                 end
            %
            %             end
            if list_EdgeFits
                disp(sprintf('mean edge slope of %5.0f good fits = %1.2e +/- %1.2e... mean OutFlagVal = %1.2e +/- %1.2e', length(Sgv), mean(Sgv), std(Sgv), mean(Sgf),std(Sgf)))
                disp(sprintf('mode edge slope of %5.0f good fits = %1.2e ... mode OutFlagVal = %1.2e', length(Sgv), mode(Sgv), mode(Sgf)))
                disp('-------------------')
                disp(sprintf('mean edge slope of %5.0f bad fits = %1.2e +/- %1.2e... mean OutFlagVal = %1.2e +/- %1.2e', length(Sbv), mean(Sbv), std(Sbv), mean(Sbf),std(Sbf)))
                disp(sprintf('mode edge slope of %5.0f bad fits = %1.2e ... mode OutFlagVal = %1.2e', length(Sbv), mode(Sbv), mode(Sbf)))
                disp('-------------------')
                disp(sprintf('mean edge slope of all data = %1.2e +/- %1.2e... mean OutFlagVal = %1.2e +/- %1.2e', mean(Out_slope), std(Out_slope), mean(Out_flagVAL(:,2)), std(Out_flagVAL(:,2))))
                disp(sprintf('mode edge slope of all data = %1.2e... mean OutFlagVal = %1.2e', mode(Out_slope), mode(Out_flagVAL(:,2))))
            end
            
            if plot_CoordXfer
                OutXfer2 = OutXfer';
                %figure(13)
                figure(48)
                subplot(6,1,1)
                imagesc(OutXfer2);
                colormap(gray)
                axis('off')
                %axis('image equal')
                xlabel('Angle from zero to 360 degrees')
                ylabel('Radius from lesion center (top) to periphery (bottom)')
                title('Angular brightness map' )
                
                %subplot(1,2,2)
                % imagesc(dat4SHOW(:,:,3))
                %colormap(gray)
                %axis('off')
                %axis('image equal')
                %title('Original Image' )
                
            end
            
            Var_Rtheta = std(Rtheta)/mean(Rtheta);  % METRIC
            Rtheta_Diff = Rtheta; % init
            
            for i_checkit = 1:length(Rtheta_Diff)-1
                Rtheta_Diff(i_checkit) = Rtheta_Diff(i_checkit) - Rtheta_Diff(i_checkit+1);
            end
            
            Rtheta_Diff (length(Rtheta_Diff)) = Rtheta_Diff (length(Rtheta_Diff)-1);
            Rtheta_Diff = abs(Rtheta_Diff);
            sumRDiff = sum(Rtheta_Diff./Rtheta); % METRIC
            
            if i_layer == 1
                Rtheta_SPEC = zeros(length(Rtheta), Nlayers);
            end
            
            Rtheta_SPEC(:,i_layer) = Rtheta;
            
            if i_layer == Nlayers
                
                for i_checkit = 1:length(Rtheta)
                    var_RAD(i_checkit) = std(Rtheta_SPEC(i_checkit, :))/mean(Rtheta_SPEC(i_checkit, :));
                end
                
                mean_var_RAD = mean(var_RAD);  % METRIC
                std_var_RAD = std(var_RAD);  % METRIC
                
            end
            
            range_mean = (max(RADmean) - min(RADmean))/mean(RADmean);
            std_mean = std(RADmean)/mean(RADmean);
            range_std = mean(RADstd./RADmean);
            std_std = std(RADstd);
            dth = 360/length(RADmean);
            theta_plot = (1:length(RADmean))*dth;
            Steve_U = (2:length(theta_plot)-2);
            
            figure(1)
            subplot(3,3,6)
            plot(theta_plot,RADmean,[clr(i_layer) '-'],'linewidth',4)
            
            %theta_I_save(:,i_layer) = RADmean;
            %theta_st_save(:,i_layer) = RADstd;
            
            if i_layer == 1
                runningMAX = max(RADmean);
            else 
                try
                    runningMAX = max([runningMAX max(RADmean)]); % for case where red shifted image yielded results
                catch
                    runningMAX = max(RADmean);
                end
            end
            
            hold('on')
            text(0-5,0,'0')
            text(90-10,0,'\pi/2')
            text(180-5,0,'\pi')
            text(270-15,0,'3\pi/2')
            text(360-10,0,'2\pi')
            
            axis('off')
            title('Angular Brightness','fontsize',16)
            hold('on')
            SmoothRad = smooth(RADmean,8);
            
            for isamp = 1:length(RADmean)
                plot([isamp isamp]*dth,[SmoothRad(isamp)-RADstd(isamp) SmoothRad(isamp)+RADstd(isamp)],'k-')
                hold('on')
            end
            axis([0 360 0 runningMAX*1.2])
            ylabel('brightness')
            plot(theta_plot,SmoothRad,[clr(i_layer) '-'],'linewidth',2)
            
            %% calculae the first order derivitave numerically
            RADdir = zeros( length(theta_plot),1 );
            for isamp = 1:length(RADmean)-1
                RADdir(isamp) = abs(SmoothRad(isamp)-SmoothRad(isamp+1));
            end
            RADdir(length(SmoothRad)) = abs(SmoothRad(length(SmoothRad))-SmoothRad(1)); % Loop around!
            
            %% Condition the function that specifies the shoulder edge sharpness
            Out_slope = Out_slope/sqrt(stats.Area); % this normalizes the edge thickness to the lesion size
            RADlengths(:,i_layer) = RADlength;
            
            %axis([0 360 0 100])
            axis('off')
            nametemp = name;
            JustNum = str2double(nametemp(2:5));
            
            %         Matrix_Out(i_layer,cnt+i,1) = range_mean;
            Matrix_Out_tmp = zeros(Nlayers,1,Nvars);
            
            Matrix_Out_tmp(i_layer,1,1) = range_mean;  % 1 used F(lambda) %done
            Matrix_Out_tmp(i_layer,1,2) = std_mean;  % 2 used  F(lambda) %done
            Matrix_Out_tmp(i_layer,1,3) = range_std;  % 3 used  F(lambda) %done
            Matrix_Out_tmp(i_layer,1,4) = std_std;  % used 4  F(lambda) %done
            
            Matrix_Out_tmp(i_layer,1,6) = std(RADdir(Steve_U));  % 6 used  F(lambda) %done
            Matrix_Out_tmp(i_layer,1,7) = sum(RADdir(Steve_U));  % 7 used  sum of hotspots (change(angle)) %done
            Matrix_Out_tmp(i_layer,1,8) = Big_Asym/stats.Eccentricity;  %branches %NOT USED F(lambda) %done
            Matrix_Out_tmp(i_layer,1,9) = Asym(i);  % branches NOT USED  F(lambda) %done
            Matrix_Out_tmp(i_layer,1,10) = brd;  %  8 used F(lambda) %done
            Matrix_Out_tmp(i_layer,1,11) = mean(Out_slope); % 9 used F(lambda) %done
            Matrix_Out_tmp(i_layer,1,12) = std(Out_slope)/mean(Out_slope);  %done % Changed for maybe (NOT USED F(lambda)
            Matrix_Out_tmp(i_layer,1,13) = fracDIM;  % 10 used F(lambda) %done
            Matrix_Out_tmp(i_layer,1,14) = colorVAR;  % 11 used F(lambda) %done
            Matrix_Out_tmp(i_layer,1,15) = mean(Out_flagVAL(:,2));  % 12 used F(lambda) %done
            Matrix_Out_tmp(i_layer,1,16) = std(Out_flagVAL(:,2));    % 13 used %done
            
            Matrix_Out_tmp(i_layer,1,19) = varPOP_length; % branches 14 used   CHANGED
            Matrix_Out_tmp(i_layer,1,20) = varMEAN_bright; % branches 15 used  CHANGED
            Matrix_Out_tmp(i_layer,1,21) = varOVERALL_bright; % branches 16 used  CHANGED
            
            Matrix_Out_tmp(i_layer,1,24) = Connectedness;  % branches 17 used
            
%             Matrix_Out_tmp(i_layer,1,25) = 0; % 18 used
%             Matrix_Out_tmp(i_layer,1,26) = 0; % NOT USED
%             Matrix_Out_tmp(i_layer,1,27) = 0; % 19 used
            
            Matrix_Out_tmp(i_layer,1,25) = rangeLen; % branches 18 used
            Matrix_Out_tmp(i_layer,1,26) = rangeClr; % branches NOT USED
            Matrix_Out_tmp(i_layer,1,27) = rangeNum; % branches 19 used
            Matrix_Out_tmp(i_layer,1,28) = sum(sum(mole)); % 20 used
            Matrix_Out_tmp(i_layer,1,29) = Var_Rtheta; % done
            Matrix_Out_tmp(i_layer,1,30) = sumRDiff; %done
            Matrix_Out_tmp(i_layer,1,31) = fracDIM_PN; % done
            Matrix_Out_tmp(i_layer,1,32) = length(Sbv)/(length(Sbv)+length(Sgv)); % done New, relatve number bad edge fits
            
            Matrix_Out_tmp(i_layer,1,33) = mean(Sgv); % done New mean edge slope of good fits
            Matrix_Out_tmp(i_layer,1,34) = std(Sgv); % done New std edge slope of good fits
            Matrix_Out_tmp(i_layer,1,35) = mean(Sgf); % done New mean edge-fit fVal of good fits
            Matrix_Out_tmp(i_layer,1,36) = std(Sgf); % done New std dge-fit fVal of good fits
            Matrix_Out_tmp(i_layer,1,37) = mode(Sgv); % done New mode edge slope of good fits
            Matrix_Out_tmp(i_layer,1,38) = mode(Sgf); % done New mode edge slope of good fits
            % Takeout Matrix_Out_tmp(i_layer,1,39) = mean(Sbv); % New mean edge slope of bad fits
            % Takeout Matrix_Out_tmp(i_layer,1,40) = std(Sbv); % New std edge slope of bad fits
            
            %             if i == 104
            %                 disp(sprintf('Coulmn 40, std(Sbv) = %0.5f',layerData(Sbv)))
            %                 disp(sprintf('i_layer = %3.0f',i_layer))
            %                 %kjh
            %             end
            
            % Takeout Matrix_Out_tmp(i_layer,1,41) = mean(Sbf); % New mean edge-fit fVal of bad fits
            % Takeout Matrix_Out_tmp(i_layer,1,42) = std(Sbf); % New std dge-fit fVal of bad fits
            % Takeout Matrix_Out_tmp(i_layer,1,43) = mode(Sbv); % New mode edge slope of bad fits
            % Takeout Matrix_Out_tmp(i_layer,1,44) = mode(Sbf); % New mode edge slope of bad fits
            Matrix_Out_tmp(i_layer,1,45) = mode(Out_slope); % done New mode to go with metric 16,17
            Matrix_Out_tmp(i_layer,1,46) = mode(Out_flagVAL(:,2))  ; % done New mode to go with metric 16,17
            Matrix_Out_tmp(i_layer,1,47) = flag2seg(i_layer); %done % New mode to go with metric 16,17
            
            % Takeout Matrix_Out_tmp(i_layer,1,48) = difcrestn;
            % Takeout Matrix_Out_tmp(i_layer,1,49) = maxdifcrestn;
            % Takeout Matrix_Out_tmp(i_layer,1,50) = mincrestloc;
            % Takeout Matrix_Out_tmp(i_layer,1,51) = maxabsmiddifcrestloc;
            Matrix_Out_tmp(i_layer,1,52) = numcrest; % done
            % Takeout Matrix_Out_tmp(i_layer,1,53) = mcrestloc;
            % Takeout Matrix_Out_tmp(i_layer,1,54) = mdifcrestloc;
            % Takeout Matrix_Out_tmp(i_layer,1,55) = absmiddifcrestloc;
            % Takeout Matrix_Out_tmp(i_layer,1,56) = eccentricity;
            % Takeout Matrix_Out_tmp(i_layer,1,57) = ellerr;
            %             Matrix_Out_tmp(i_layer,1,58) = mean(Sgvlinmid); % New mean edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,59) = std(Sgvlinmid); % New std edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,60) = mean(Sgflinmid); % New mean edge-fit fVal of good fits
            %             Matrix_Out_tmp(i_layer,1,61) = std(Sgflinmid); % New std dge-fit fVal of good fits
            %             Matrix_Out_tmp(i_layer,1,62) = mode(Sgvlinmid); % New mode edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,63) = mode(Sgflinmid); % New mode edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,64) = mean(Sgvlinbeg); % New mean edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,65) = std(Sgvlinbeg); % New std edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,66) = mean(Sgflinbeg); % New mean edge-fit fVal of good fits
            %             Matrix_Out_tmp(i_layer,1,67) = std(Sgflinbeg); % New std dge-fit fVal of good fits
            %             Matrix_Out_tmp(i_layer,1,68) = mode(Sgvlinbeg); % New mode edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,69) = mode(Sgflinbeg); % New mode edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,70) = mean(Sgvlinend); % New mean edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,71) = std(Sgvlinend); % New std edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,72) = mean(Sgflinend); % New mean edge-fit fVal of good fits
            %             Matrix_Out_tmp(i_layer,1,73) = std(Sgflinend); % New std dge-fit fVal of good fits
            %             Matrix_Out_tmp(i_layer,1,74) = mode(Sgvlinend); % New mode edge slope of good fits
            %             Matrix_Out_tmp(i_layer,1,75) = mode(Sgflinend); % New mode edge slope of good fits
            %
            if i_layer == Nlayers
                %% make color image
                sliceR = dat4MATH(:,:,1);
                sliceG = dat4MATH(:,:,2);
                sliceB = dat4MATH(:,:,3);
                
                % LukUp = [ % from Mete et al. CMIG-1138 2012
                % %   minR    minG    minB    maxR    maxG    maxB
                %     222     208     195     255     255     255
                %     225     11      0       243     25      7
                %     149     81      16      186     108     42
                %     69      11      10      101     30      8
                %     69      79      68      92      110     120
                %     0       0       0       29      26      11 ];
                
                % Training set
                % 87 is D07 small blue/grey in center
                % 92 is D12 large blue/grey in center
                % 125 is D45 large blue/grey in center
                % 137 is D57 small blue/gray spot lower left
                % 145 is D65 small blue/gray EVERYWHERE
                % 158 is D91 small blue/gray EVERYWHERE
                % indyE = [87 92 125 137 145 158]
                
                Colors = [ % wikipedia
                    40      26      13      % dark-brown
                    71      40      11      % light-brown
                    0       0       0       % black
                    100     0       0       % red
                    40      60      80      % blue-gray
                    100     100     100];   % white
                
                Colors = Colors/100*256;
                
                Rm = [ % ratio mean   [r/b r/g b/g]  mean_rb_rg_bg
                    1.9006    2.0193    1.0656      % dark-brown
                    1.7247    1.6208    0.9431      % light brown
                    
                    0.4648    0.7536    1.7404      % black
                    1.8058    1.9820    1.1040      % red
                    %0.8286    1.0834    1.3094      % blue-gray
                    1.2598    1.3210    1.0515
                    %1.2649    1.2254    0.9768      % blue-gray
                    0.9243    1.2008    1.2998      % white
                    ];
                
                Rs = 3*[ % std
                    0.1429    0.1344    0.0721      % dark brown
                    0.1521    0.0877    0.0479      % light brown
                    
                    0.1841    0.2127    0.3964      % black
                    0.2301    0.2032    0.0939      % red
                    %0.0702    0.0780    0.0464      % blue-gray
                    0.1143    0.0829    0.0436
                    %0.0666    0.0440    0.0293      % blue-gray
                    0.0342    0.0294    0.0257      % white
                    ];
                
                Rs(1,:) = Rs(1,:)*8;
                Rs(2,:) = Rs(2,:)*3;
                Rs(4,:) = Rs(4,:)/2;
                %Rs(5,:) = Rs(5,:);
                Rs(6,:) = Rs(6,:)*2;
                
                COLORimg = zeros(Ny,Nx,3)+256;  % make solid color "painting"
                ColorMaps = zeros(Ny,Nx,length(Rs));
                FLAGS = zeros(6,1);
                BlueFlag = -1;
                
                for indexY = 1:Ny           % scan image in y
                    for indexX = 1:Nx       % scan image in x
                        r_b = sliceR(indexY,indexX)/sliceB(indexY,indexX); % ratio of red to blue
                        r_g = sliceR(indexY,indexX)/sliceG(indexY,indexX);
                        b_g = sliceB(indexY,indexX)/sliceG(indexY,indexX);
                        for indexZ = 1:length(Rs)  % test to see if current pixel is each of 6 colors
                            % indexing: Rm(indexZ,2)  indexZ->row (putatuve color)
                            % 2-> pixeld ratio (r/g)
                            if      r_g <= Rm(indexZ,2)+Rs(indexZ,2)  && r_g >= Rm(indexZ,2)-Rs(indexZ,2) ...
                                    && r_b <= Rm(indexZ,1)+Rs(indexZ,1)  && r_b >= Rm(indexZ,1)-Rs(indexZ,1) ...
                                    && b_g <= Rm(indexZ,3)+Rs(indexZ,3)  && b_g >= Rm(indexZ,3)-Rs(indexZ,3)
                                if mole(indexY,indexX) % if pixel is inside lesion
                                    ColorMaps(indexY,indexX,indexZ) = 1;
                                    COLORimg(indexY,indexX,1) = Colors(indexZ,1);
                                    COLORimg(indexY,indexX,2) = Colors(indexZ,2);
                                    COLORimg(indexY,indexX,3) = Colors(indexZ,3);
                                    FLAGS(indexZ) = 1;
                                end
                            end
                        end
                    end
                end
                
                if sum(sum(ColorMaps(:,:,5))) > sum(sum(mole))/1000 | sum(sum(ColorMaps(:,:,6))) > sum(sum(mole))/1000
                    BlueFlag = 1;  % METRIC
                else
                    BlueFlag = 0;
                end
                
                figure(6);
                subplot(1,3,3)
                imagesc(COLORimg/max(max(max(COLORimg))));
                axis equal image
                axis('off')
                
                i_colors = 0;
                sums_color = sum(sum(ColorMaps));
                for i_check_colors = 1:length(Colors)
                    if sums_color(i_check_colors)
                        i_colors = i_colors + 1;
                    end
                end
                title(sprintf('%1.0f colors found',i_colors));
                
                %%
                % Map to CIE LAB color space
                mblack = mean(blackwb0);
                stdblack = max(blackwb0) - min(blackwb0);
                mbluegray = mean(bluegraywb0);
                stdbluegray = max(bluegraywb0) - min(bluegraywb0);
                mdbrown = mean(dbrownwb0);
                stddbrown = max(dbrownwb0) - min(dbrownwb0);
                mlbrown = mean(lbrownwb0);
                stdlbrown = max(lbrownwb0) - min(lbrownwb0);
                mwhite = mean(whitewb0);
                stdwhite = max(whitewb0) - min(whitewb0);
                mred = mean(redwb0);
                stdred = max(redwb0) - min(redwb0);
                
                mcolors = [mdbrown; mlbrown; mblack; mred; mbluegray; mwhite];
                stdcolors = [stddbrown; stdlbrown; stdblack; stdred; stdbluegray; stdwhite];
                
                ColorsCIELAB = rgb2lab(mcolors);
                temp = mcolors + stdcolors;
                stdCIELAB = sum((rgb2lab(temp) - ColorsCIELAB).^2,2);
                
                BlueFlagCIELAB = -1;
                mole3 = repmat(mole, [1 1 3]);
                moleCIE = datt./255.*mole3;
                xyz = rgb2lab(moleCIE);
                xyz2 = reshape(xyz, [size(xyz,1).*size(xyz,2), 3]);
                coldiff = zeros(size(xyz2,1),size(ColorsCIELAB,1));
                for melcol = 1:size(ColorsCIELAB,1)
                    coldiff(:,melcol) = sum(gsubtract(xyz2, ColorsCIELAB(melcol,:)).^2,2);
                end
                [diff, predcol] = min(coldiff,[],2);
                CIELABmap = zeros(Ny, Nx, 3);
                for indexX = 1:Nx
                    for indexY = 1:Ny
                        CIELABmap(indexY,indexX,1) = Colors(predcol((indexX-1).*size(CIELABmap,1)+indexY),1);
                        CIELABmap(indexY,indexX,2) = Colors(predcol((indexX-1).*size(CIELABmap,1)+indexY),2);
                        CIELABmap(indexY,indexX,3) = Colors(predcol((indexX-1).*size(CIELABmap,1)+indexY),3);
                    end
                end
                
                %     figure;
                %     imagesc(CIELABmap/max(max(max(CIELABmap))));
                %     axis equal image
                %     axis('off')
                
                CIELABmap = zeros(Ny, Nx, 3);
                stdCIELAB = stdCIELAB';
                stdCIELAB = repmat(stdCIELAB, [size(coldiff,1), 1]);
                temp = coldiff <= stdCIELAB;
                [x, predcol] = find(temp == 1);
                tempX = ceil(x./Ny);
                tempY = mod(x,Ny);
                tempY(tempY == 0) = Ny;
                for indexX = 1:length(tempX)
                    temp1 = tempY(indexX);
                    temp2 = tempX(indexX);
                    CIELABmap(temp1,temp2,1) = Colors(predcol(indexX),1);
                    CIELABmap(temp1,temp2,2) = Colors(predcol(indexX),2);
                    CIELABmap(temp1,temp2,3) = Colors(predcol(indexX),3);
                end
                
                figure(6)
                subplot(1,3,2)
                imagesc(CIELABmap/max(max(max(CIELABmap))));
                title(sprintf('Amanda:$1.0f colors found',3))
                axis equal image
                axis('off')
                
                i_colorsCIELAB = length(unique(predcol));
                if ~isempty(find(predcol == 5))
                    BlueFlagCIELAB = 1;
                else
                    BlueFlagCIELAB = 0;
                end
                
                Matrix_Out_tmp(i_layer,1,78) = clroff; %done adjusted maybe (NOT USED
                Matrix_Out_tmp(i_layer,1,79) =sz_dif; %done 21 used
                Matrix_Out_tmp(i_layer,1,80) =sz_dif_max; %done 22 used
                Matrix_Out_tmp(i_layer,1,81) =mean_var_RAD; % done
                Matrix_Out_tmp(i_layer,1,82) =std_var_RAD; % done
                %Matrix_Out_tmp(i_layer,1,83) = i_colors; % brand new, [0-6]
                %INTAGER METRIC
                Matrix_Out_tmp(i_layer,1,84) = BlueFlag; % done brand new, BINARY METRIC  REPLACED WITH i_colorsCIELAB (86)
                Matrix_Out_tmp(i_layer,1,85) = sum(flag2seg); % done brand new, BINARY METRIC
                Matrix_Out_tmp(i_layer,1,86) = i_colorsCIELAB; %done
                Matrix_Out_tmp(i_layer,1,87) =BlueFlagCIELAB; %done
                
                %Find_Vessles
                figure(1)
                dummyIM = ones(150,1000,3);
                subplot(3,3,7)
                imagesc(dummyIM)
                axis('off')
                for i_SR = 1:round(Nvars/4)
                    text(1,15*(i_SR-1)-40,sprintf('V%2.0d %1.2e',i_SR, Matrix_Out_tmp(i_layer,1,i_SR) ))
                end
                
                for i_SR = round(Nvars/4)+1:Nvars/2
                    text(600,15*(i_SR-1 - round(Nvars/4))-40,sprintf('V%2.0d= %1.2e',i_SR, Matrix_Out_tmp(i_layer,1,i_SR) ))
                end
                
                subplot(3,3,8)
                imagesc(dummyIM)
                axis('off')
                for i_SR = round(Nvars/2)+1:round(3*Nvars/4)
                    text(1,15*(i_SR-1-round(Nvars/2))-40,sprintf('V%2.0d %1.2e',i_SR, Matrix_Out_tmp(i_layer,1,i_SR) ))
                end
                
                for i_SR = round(3*Nvars/4)+1:Nvars
                    text(600,15*(i_SR-1 - round(3*Nvars/4))-40,sprintf('V%2.0d= %1.2e',i_SR, Matrix_Out_tmp(i_layer,1,i_SR) ))
                end
                
                figure(1)
                subplot(3,3,9)
                tempTHplot = (1:length(Assym2))./length(Assym2)*360;
                plot(tempTHplot,Assym2*50,'kx','markersize', 8)
                hold('on')
                plot(tempTHplot(i_Asym),Assym2(i_Asym)*50,'ro','markerfacecolor','r')
                plot(tempTHplot(i_sym),Assym2(i_sym)*50,'bo','markerfacecolor','b')
                plot([1:length(RADlengths)]/length(RADlengths)*max(tempTHplot),RADlengths(:,1)/max(max(RADlengths))*max(Assym2*50)*2,'r-')
                plot([1:length(RADlengths)]/length(RADlengths)*max(tempTHplot),RADlengths(:,2)/max(max(RADlengths))*max(Assym2*50)*2,'g-')
                plot([1:length(RADlengths)]/length(RADlengths)*max(tempTHplot),RADlengths(:,3)/max(max(RADlengths))*max(Assym2*50)*2,'b-')
                axis([0 360 0 max(Assym2*50)*2])
                axis('off')
                text(0-5,-5,'0')
                text(90-10,-5,'\pi/2')
                text(180-5,-5,'\pi')
                text(270-15,-5,'3\pi/2')
                text(360-10,-5,'2\pi')
                
                if plotON_Topo
                    H = fspecial('average', [10 10]);
                    Topo_sm = imfilter(-Topo, H);
                    X_surf = 1:Nx;
                    Y_surf = 1:Ny;
                    Z_surf = -1:min(min(Topo));
                    figure(10);clf();
                    surf(Topo_sm.*(Topo>0))
                    hold on
                   % surface( imagesc(dat4SHOW/256), [0 1])

                    %                 pause(1)
                    %                 direction = [1 0 0];
                    %                 rotate(gcf,direction,10)
                    %                 drawnow
                    
                    % h = rotate3d;
                    % setAllowAxesRotate(h,gca,false);
                    
%                     for iii = 1:360
%                         view(iii,45)
%                         pause(0.05)
%                     end
%                     
                end
                
                if BrnchFrqAnalys
                    freq_comps(i).img = freqMAP_out;
                end
                
            end
            
            
            %             if i_layer == 1
            %                 Matrix_Out1(:,i+cnt,:) = Matrix_Out_tmp;
            %             end
            %             if i_layer == 2
            %                 Matrix_Out2(:,i+cnt,:) = Matrix_Out_tmp;
            %             end
            %             if i_layer == 3
            %                 Matrix_Out3(:,i+cnt,:) = Matrix_Out_tmp;
            %             end
            %
            if i_layer == 1
                Matrix_Out1(:,i,:) = Matrix_Out_tmp;
            end
            if i_layer == 2
                Matrix_Out2(:,i,:) = Matrix_Out_tmp;
            end
            if i_layer == 3
                Matrix_Out3(:,i,:) = Matrix_Out_tmp;
            end
            
            figure(1)
            subplot(3,3,4)
            lesion_ctr = round(stats.Centroid);
            plot(lesion_ctr(1),lesion_ctr(2),[clr(i_layer) '*'])  %Lookhere
            
            subplot(3,3,5)
            imagesc(datt/256)
            hold('on')
            axis equal image
            axis('off')
            
            Diagnoses(i,10) = 0; % record one or more channel segmentation su
            catch
                Diagnoses(i,10) = 1; % record one or more channel segmentation errors
                 
            end
            
        end
        
        figure(48)
        subplot(6,1,2)
        % for i_layer = 1:Nlayers
        %         plot(theta_plot',theta_I_save(:,i_layer),[clr(i_layer) '-'],'linewidth',2)
        %         hold on
        %
        %         for i_doit = 1:length(theta_plot)
        %             plot([theta_plot(i_doit)  theta_plot(i_doit)],[    (theta_I_save(i_doit,i_layer) -theta_st_save(i_doit,i_layer))      (theta_I_save(i_doit,i_layer) +theta_st_save(i_doit,i_layer))       ],'k-')
        %         end
        %     end
        %     axis([0 360 0 max(max(theta_I_save))])
        
        %          std_std = std(RADstd);
        %         dth = 360/length(RADmean);
        %         theta_plot = (1:length(RADmean))*dth;
        %         Steve_U = (2:length(theta_plot)-2);
        %
        %         figure(1)
        %         subplot(3,3,6)
        %         plot(theta_plot,RADmean,[clr(i_layer) '-'],'linewidth',4)
        %
        %         theta_I_save(:,i) = RADmean;
        %         theta_st_save(:,i) = RADstd;
        
        
        
        
        
        % grant: how do I get the figures saved?
        
        if SaveOutFigs
            cd(targetDIR)
            if ~isdir('ResultsImages');
                mkdir ResultsImages;
            end
            cd ResultsImages
            figs2PRINT = [1 ]
            for i_print_it = 1:length(figs2PRINT)
                fig_prnt_no = figs2PRINT(i_print_it);
                name = sprintf(['Out', num2str(fig_prnt_no),name]);
                set(figure(fig_prnt_no),'position',[200  100  900  550])
                print(figure(figs2PRINT(i_print_it)),'-djpeg','-r600',name);
            end
            cd(homeDIR)
        end
        
        catch
           errim(i,1) = -1;
           disp(sprintf('ERROR, skipping lesion %3.0f', i));
        end
        
            %else
             %  disp(sprintf('Skipped Lesion %4.0 because it was not a melanoma or nevus',i))
            %end % if melanoma or nevus
    disp(sprintf('Finished Lesion %3.0f of %3.0f in %4.0f[s]', i,Nims,toc));
end

Matrix_Out_Sum = Matrix_Out1 + Matrix_Out2 + Matrix_Out3;
[NzM NyM NxM] = size(Matrix_Out_Sum);
for iz = 1:NzM
    for iy = 1:NyM
        for ix = 1:NxM
            % Dan removed the following if statement because there were no 
            % BlueFlagCIELAB = 1 signals in a set of nevi hence Matrix_Out
            % wasn't being allocated to the full size (BlueFlagCIELAB was last
            % coulmn)
            
            %if abs(Matrix_Out_Sum(iz,iy,ix)) > 0;  
                Matrix_Out(iz,iy,ix) = Matrix_Out_Sum(iz,iy,ix);
            %end
        end
    end
end

try
    Diagnoses(errim == -1, 9) = 1;
    Diagnoses(errim == -1, 20) = 1;
catch
end

if using_parfor
    cd ..
    cd Data
    cd(folderDIR)
else
    cd(DataFolder)
end
save GSdiag Diagnoses
save Matrix_mIBs Matrix_Out

if using_parfor
    delete(gcp);
end
Matrix_Write = zeros(Nims,Nvars);

toc

for i_layer = 1:Nlayers
    Matrix_Write = shiftdim(Matrix_Out(i_layer,:,:),1);
    xlswrite(sprintf(['AutoSaveResult' num2str(i_layer)]),Matrix_Write)
end


if ~isdir('Track_mIBs')
    mkdir('Track_mIBs')
end

cd Track_mIBs
for i_layer = 1:Nlayers
    Matrix_Write = shiftdim(Matrix_Out(i_layer,:,:),1);
    xlswrite(sprintf([date '_AutoSaveResult' num2str(i_layer)]),Matrix_Write)
end
statement = ['save ' date '_GSdiag Diagnoses'];
eval(statement)
statement = ['save ' date '_Matrix_mIBs Matrix_Out'];
eval(statement)

cd ..
beep; pause(0.2);beep; pause(0.2);beep; pause(0.2);