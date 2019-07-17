function [Out_img] = GausTail(In_img, PerimPts)


Out_img = In_img;

[Ny Nx] = size(In_img);
CCC = Nx/10;
for i_x  = 1:Nx
    for i_y = 1:Ny
        if In_img(i_y,i_x) == 0
            dist = Nx;
            for i_seek = 1:length(PerimPts)
                dist2 = sqrt(  (PerimPts(1,i_seek)-i_x)^2 + (PerimPts(2,i_seek)-i_y)^2 );
                if dist2 < dist
                    dist = dist2;
                    CloseX = PerimPts(1);
                    CloseY = PerimPts(2);
                end
            end
            Out_img(i_y,i_x) = In_img(CloseY,CloseX) * exp(-dist^2/CCC^2);
        end
    end
end
% figure(76)
% subplot(2,1,1)
% imagesc(In_img)
% colormap gray
% colorbar
% 
% subplot(2,1,2)
% imagesc(Out_img)
% colormap gray
% colorbar

end