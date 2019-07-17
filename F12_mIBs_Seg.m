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
end

mole = imfill(mole,'holes');
%stop
if plotON_Thresh
    subplot(3,3,4)
    imagesc(mole)
    axis('image')
    title('imfill')
    axis('off')
    colormap(gray)
end

mole = imerode(mole,seD);

if plotON_Thresh
    subplot(3,3,5)
    imagesc(mole)
    axis('image')
    title('imerode')
    axis('off')
    colormap(gray)
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