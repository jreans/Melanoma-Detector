function [skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts] = TraceBranch(skel2,lastpix,Current_branch_length,BrPxLst,branchPts,endPts)
% TraceBranch.m
jam_branch = 1;
while jam_branch
    nextpix = [0 0];
    skel2(lastpix(1),lastpix(2)) = 0; 
    if  skel2(lastpix(1)-1,lastpix(2)) 
        nextpix = [lastpix(1)-1 lastpix(2)];
    end
    if  skel2(lastpix(1),lastpix(2)+1) 
        nextpix = [lastpix(1) lastpix(2)+1];
    end
    if  skel2(lastpix(1)+1,lastpix(2))
         nextpix = [lastpix(1)+1 lastpix(2)];
    end
    if  skel2(lastpix(1),lastpix(2)-1)
        nextpix = [lastpix(1) lastpix(2)-1];
    end       
    if sum(nextpix) == 0
            if  skel2(lastpix(1)-1,lastpix(2)+1) 
                nextpix = [lastpix(1)-1 lastpix(2)+1];
            end
            if  skel2(lastpix(1)+1,lastpix(2)+1) 
                nextpix = [lastpix(1)+1 lastpix(2)+1];
            end
            if  skel2(lastpix(1)+1,lastpix(2)-1)
                 nextpix = [lastpix(1)+1 lastpix(2)-1];
            end
            if  skel2(lastpix(1)-1,lastpix(2)-1)
                nextpix = [lastpix(1)-1 lastpix(2)-1];
            end    
    end
    Current_branch_length = Current_branch_length + 1;
    BrPxLst(Current_branch_length,1) = nextpix(1);
    BrPxLst(Current_branch_length,2) = nextpix(2);
    for i_look2 = 1:length(branchPts)  % reach an end of branch???
        if nextpix == branchPts(i_look2,:)
            jam_branch = 0;
        end
    end
    for i_look2 = 1:length(endPts)  % reach an end of branch???
        if nextpix == endPts(i_look2,:) 
            jam_branch = 0;
        end
    end
    if sum(nextpix) == 0 %this is the case for branch of length 1
        jam_branch = 0;
    end
    figure(88)
    subplot(3,6,17)
    plot(lastpix(2),lastpix(1),'y.')
    lastpix = nextpix;
end
