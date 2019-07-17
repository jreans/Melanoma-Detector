% Run4_Pvalues.m

clear all
close all

homedir = pwd;
load CurrentWorkingDir
cd Data
cd(folderDIR)
load Matrix_mIBs
load GSdiag % overwrites Diagnoses of Matrix_OutProcessed

SaveName = 'Mel_vs_Nevus';
if SaveName == 'Mel_vs_Nevus'
    LookType1 = 0;
    LookType2 = 1;
end

if strcmp(folderDIR,'S2_Puig')   % special case for puig, grabs the correct diagnoses
    Diagnoses(:,2) = Diagnoses(:,15);
end
if strcmp(folderDIR,'S11_UCI_OHSU')   % special case for puig, grabs the correct diagnoses
    Diagnoses(:,2) = Diagnoses(:,1);
end

save GSdiag Diagnoses

[folderDIR sprintf(' %3.0f images worked of %3.0f ', length(Diagnoses) - sum(Diagnoses(:,20)), length(Diagnoses) )]
[Nc Nl Nb] = size(Matrix_Out); % colors lesions biomarkers

Diagnoses(:,20) = 0; % over-ride original [bad flags] and set by imaging biomarker availability below
Diagnoses(:,19) = 0;

RGBname = 'RGB';

i_nothing = 0; % INITIALIZATIONS
i_NotAll = 0;
i_isNaN = 0;
n_nothing = 0;
n_NotAll = 0;
n_NaN = 0;

for i_l = 1:Nl
    Net_Biom = sum(shiftdim(Matrix_Out(:,i_l,:),2)); % yields the sum of the biomarker value for the red, green, blue channels
    if min(abs(Net_Biom)) > 0  % has at least 1 mIB value in at least one color channel
        for i_b = 1:Nb
            temp_look = abs(Matrix_Out(:,i_l,i_b));
        end
        Diagnoses(i_l,19) = 0;
    elseif sum(Net_Biom) == 0
        i_nothing = i_nothing + 1;
        n_nothing(i_nothing) = i_l;
        disp(sprintf('Lesion %4.0f had no biomarkers in any channel', i_l))
        Diagnoses(i_l,19) = 1;
    elseif sum((Net_Biom==0)) < 3 & sum((Net_Biom==0)) >0
        i_NotAll= i_NotAll + 1;
        n_NotAll(i_NotAll) = i_l;
        disp([sprintf('Lesion %4.0f had no biomarkers in the ', i_l) RGBname(Net_Biom==0) ' channel(s)'])
        Diagnoses(i_l,19) = 1;
    else
        i_isNaN = i_isNaN + 1;
        n_NaN(i_isNaN) = i_l;
        disp([sprintf('Lesion %4.0f had NaNs in the ', i_l) RGBname(isnan(Net_Biom)) ' channel(s)'])
        Diagnoses(i_l,19) = 1;
    end
end



disp('--------------------------------------------')
% disp(['Lesions with no mIB vals: ' num2str(n_nothing)])
% disp(['Lesions with not all mIB vals: ' num2str(n_NotAll)])
% disp(['Lesions with NaN mIB vals: ' num2str(n_NaN)])

disp([sprintf('%4.0f Lesions with no mIB vals: ', length(n_nothing) ) num2str(n_nothing) ])
disp([sprintf('%4.0f Lesions with not all mIB vals: ', length(n_NotAll) ) num2str(n_NotAll) ])
disp([sprintf('%4.0f Lesions with NaN mIB vals: ', length(n_NaN) ) num2str(n_NaN) ])

Matrix_TmpR = shiftdim(Matrix_Out(1,:,:),1);
Matrix_TmpG = shiftdim(Matrix_Out(2,:,:),1);
Matrix_TmpB = shiftdim(Matrix_Out(3,:,:),1);

cd(homedir)

Matrix_Out2 = Matrix_Out;  % initialize to later strip errors out of
Diagnoses2 = Diagnoses; % initialize to later strip errors out of
diag_only = Diagnoses(:,2);

X1 = diag_only; % initialize to later strip errors out of

[blah Size_MO blah2] = size(Matrix_Out2);
Size_Diag = length(X1);

if Size_Diag > Size_MO  % ?????????????????????????????
    X1(Size_Diag) = [];
end

[xx1 Nl xx2] = size(Matrix_Out2); % colors lesions biomarkers

if strcmp(folderDIR,'S2_Puig')   % special case for puig
    cd(DataFolder)
    load ManualChecks
    for i = 1:length(ManualCheck)
        if ManualCheck(i) == 0
            ManualCheck(i) = 1;
        elseif ManualCheck(i) == 1
            ManualCheck(i) = 0;
        elseif ManualCheck(i) < 0
            ManualCheck(i) = 1;
        end
    end
    Diagnoses(:,20) = ManualCheck';
end

N_mIBs_didnt_work = i_nothing + i_NotAll + i_isNaN;
N_autocrop_didnt_work = sum(Diagnoses(:,20));

for i = 1:Nl  % remove errors  But Diagnoses(:,1) is the original image number
    i_neg = Nl-i+1;
    kil_it = 0;
    if Diagnoses(i_neg,20) | Diagnoses(i_neg,19)
        kil_it = 1;
    end
    if Diagnoses(i_neg,2) ~=LookType1 & Diagnoses(i_neg,2) ~=LookType2  % melanoma or nevus onlys
        kil_it = 1;
    end
    if kil_it == 1
        Diagnoses2(i_neg,:) = [];
        Matrix_Out2(:,i_neg,:) = [];
        X1(i_neg) = [];
    end
end

disp(sprintf('Found %3.0f of lesion type 1 and %3.0f of lesion type 2', length(X1(X1==LookType1)), length(X1(X1==LookType2))));

% Write out mIBs ffor Amanda where the bad lesions are stripped but all
% color channel mIBs are included.


%imID = Diagnoses(Diagnoses(:,20) == 0 & Diagnoses(:,19) == 0,1);
i_cnt = 0;
for i = 1:length(Diagnoses)
    if Diagnoses(i,20) == 0 & Diagnoses(i,19) == 0
        i_cnt = i_cnt+1;
        imID(i_cnt) = i;
    end
end

[Nc2 Nl2 Nb2] = size(Matrix_Out2); % colors lesions biomarkers
temp_mIBs = zeros(Nl2,1); % initialize
pval_out = zeros(Nc2, Nb2) - 1;  % initialize as -1

cd ../..

for i_c = 1:Nc
    Matrix_Tmp = shiftdim(Matrix_Out2(i_c,:,:),1);
    for i_b = 1:Nb
        clear X2
        for i_l = 1:Nl2
            X2(i_l) = Matrix_Out2(i_c,i_l,i_b);
        end
        X2 = X2';
        
        %% determine the p-value for varius input data types: binary, muliclass, and continuous
        if max(X2) == 1 & sum(round(X2) - X2) == 0 %% identify binary for exact fisher test
            i_b
            TP = 0; % Construct a confusion matrix
            TN = 0;
            FP = 0;
            FN = 0;
            for i_scan = 1:Nl2
                if X1(i_scan) == LookType2 & X2(i_scan)
                    TP = TP + 1;
                elseif ~X1(i_scan)== LookType2 & ~X2(i_scan)
                    TN = TN + 1;
                elseif ~X1(i_scan)== LookType2 & X2(i_scan)
                    FP = FP + 1;
                elseif X1(i_scan)== LookType2 & ~X2(i_scan)
                    FN = FN + 1;
                end
            end
            
            x = table([TP;FN],[FP;TN],'VariableNames',{'Mela','Nevi'},'RowNames',{'posVal','negVal'})
            [h,p,stats] = fishertest(x,'Tail','right','Alpha',0.01);
            pval_out(i_c,i_b) = p;
            
        elseif max(X2) < 10 & sum(round(X2) - X2) == 0 & abs(sum(X2)) >0 %% identify discrete catagorical
            clear results
            ALPHA = 0.05;
            X3 = X2(X1==LookType2);
            X4 = X2(X1==LookType1);
            % [p,h,stats] = ranksum(X1,X2);
            [p,h,stats] = ranksum(X3,X4);
            pval_out(i_c,i_b) = p;
            
        elseif abs(sum(round(X2) - X2)) > 0  %% identify continuous
            TST = 0; % TST - unpaired (0) or paired (1) test (default = 0).
            ALPHA = 0.05; %  ALPHA - significance level (default = 0.05).
            TAIL = 2; % TAIL - 1-tailed test (1) or 2-tailed test (2). (default = 2).
            clear results
            X3 = X2(X1==1);
            X4 = X2(X1==0);
            cd(homedir)
            results = F14_tTest(X3,X4,TST,ALPHA,TAIL);
            pval_out(i_c,i_b) = results.tpvalue;
        end
    end
end
cd(homedir)

figure(1)
if PCscreen
    set(figure(1),'position',1000*[ 0.0010    0.0383    1.2800    0.6907])
end
clr = 'rgbm';
mIBs2use = zeros(Nb2); % initialize

%%
subplot(3,1,1)

for i_c = 1:Nc2
    
    Matrix_Tmp = shiftdim(Matrix_Out2(i_c,:,:),1);
    delta = (i_c-1)*0.2;
    for i_b = 1:Nb2
        tempP = pval_out(i_c,i_b);
        if tempP > 0
            mIBs2use(i_b) = 1;
            pos = i_b+delta;
            plot([pos pos], [0 -log10(tempP)], [clr(i_c) '-'],'linewidth',2)
            hold on
            if i_c == 3
                text(i_b,-1, sprintf('%2.0f',i_b),'fontsize',6)
            end
        end
    end
end

plot([0 pos],[-log10(0.05) -log10(0.05)],'k--')
ylabel('-log10(p)')
xlabel('mIB coulmn')

mc_flag = zeros(Nb2,1); % initialize flags for multicolor metrics
for i_b = 1:Nb
    sum1 = sum(Matrix_Out2(1,:,i_b));
    sum2 = sum(Matrix_Out2(2,:,i_b));
    sum3 = sum(Matrix_Out2(3,:,i_b));
    if sum1+sum2 == 0 & sum3 >0
        mc_flag(i_b) = 1; % Identify multicolor metrics
    end
end

good = zeros(Nb2,1);  %Finds imaging biomarkers where there are all three p-values
for i_b = 1:Nb
    if pval_out(1,i_b) > 0 & pval_out(2,i_b) > 0 & pval_out(3,i_b) > 0
        good(i_b) = 1;
    end
end

i_count_mc = 0; % initialize to count multicolor metrics
for i_c = 1:Nc2
    delta = (i_c-1)*0.2;
    for i_b = 1:Nb2
        tempP = pval_out(i_c,i_b);
        if mc_flag(i_b) & i_c == Nc
            i_count_mc = i_count_mc + 1;
            MCp(1,i_count_mc) = i_b;
            MCp(2,i_count_mc) = tempP;
        end
    end
end



subplot(3,1,2)
[NPVy NPVx] = size(pval_out);
pval_out2 = zeros(NPVy+1,NPVx); % initialize to later strip out the bad mIBs
pval_out2(2:NPVy+1,:) = pval_out; % put in an extra row for mIB couplmn #
pval_out2(1,:) = [1:Nb]; % put in an extra row for mIB coulmn #



for i_b = 1:Nb  % remove errors  But pval_out2(1,:) is the original mIB coulmn
    i_neg = Nb-i_b+1;
    if min(pval_out2(:,i_neg)) <= 0
        pval_out2(:,i_neg) = [];
    end
end



[four Nb3] = size(pval_out2);

for i_c = 1:Nc2
    delta = (i_c-1)*0.2;
    for i_b = 1:Nb3
        tempP = pval_out2(i_c+1,i_b);
        plot(delta + [i_b i_b], [0 -log10(tempP)], [clr(i_c) '-'],'linewidth',3)
        hold on
        if i_c == 3
            text(i_b,-1, sprintf('%2.0f',pval_out2(1,i_b)),'fontsize',6)
        end
    end
end

for i_mc = 1:length(MCp)
    tempP = MCp(2,i_mc);
    pos2 = Nb3+i_mc;
    plot([pos2 pos2], [0 -log10(tempP)], [clr(i_c+1) '-'],'linewidth',10)
    text(pos2,-1, sprintf('%2.0f',MCp(1,i_mc)),'fontsize',6)
end

axis([0 (sum(good) + i_count_mc +1) 0 (1.2*max(max(-log10(pval_out))))])
plot([0 pos],[-log10(0.05) -log10(0.05)],'k--')
ylabel('-log10(p)')
xlabel('mIB Number')

%% Find most significant color version of each single color channel mIB
Nr = 0;
Ng = 0;
Nb = 0;
for i_b = 1:Nb3
    if pval_out2(2,i_b) == min(pval_out2(2:4,i_b))
        Nr = Nr+1;
        mIBs_Red(:,Nr) = pval_out2(:,i_b);
    end
    if pval_out2(3,i_b) == min(pval_out2(2:4,i_b))
        Ng = Ng+1;
        mIBs_Green(:,Ng) = pval_out2(:,i_b);
        %Green
    end
    if pval_out2(4,i_b) == min(pval_out2(2:4,i_b))
        Nb = Nb+1;
        mIBs_Blue(:,Nb) = pval_out2(:,i_b);
    end
end


%% create four subgroups of ordered mIBS: red, green, blue, multicolor
mIBs_Red = mIBs_Red';
try
    mIBs_Green = mIBs_Green';
catch
end
mIBs_Blue = mIBs_Blue';
MCp = MCp';

mIBs_Red_Sorted = sortrows(mIBs_Red,2);
try
    mIBs_Green_Sorted = sortrows(mIBs_Green,3);
catch
end
mIBs_Blue_Sorted = sortrows(mIBs_Blue,4);
MCp_Sorted = sortrows(MCp,2);

MCp_Sorted = MCp_Sorted';
mIBs_Red_Sorted = mIBs_Red_Sorted';
try
    mIBs_Green_Sorted = mIBs_Green_Sorted';
catch
end
mIBs_Blue_Sorted = mIBs_Blue_Sorted';


%% make connections to original paper's mIBs

LabelName(83).name = 'MC1'; % LabelName(53).name = 'MC1';
LabelName(80).name = 'MC2'; % LabelName(50).name = 'MC2';
LabelName(81).name = 'MC3'; % LabelName(51).name = 'MC3';
LabelName(84).name = 'MC4'; % LabelName(54).name = 'MC4';

LabelName(7).name = 'B1';
LabelName(4).name = 'B2';
LabelName(41).name = 'B3';
LabelName(44).name = 'B4';
LabelName(6).name = 'B5';
LabelName(15).name = 'B6';
LabelName(35).name = 'B7';
LabelName(20).name = 'B8';
LabelName(42).name = 'B9';
LabelName(2).name = 'B10';
LabelName(19).name = 'B11';
LabelName(1).name = 'B12';
LabelName(16).name = 'B13';
LabelName(36).name = 'B14';
LabelName(10).name = 'G1';
LabelName(21).name = 'B15';

LabelName(34).name = 'R1';
LabelName(13).name = 'R2';
LabelName(24).name = 'R3';
LabelName(28).name = 'R4';
LabelName(33).name = 'R5';
LabelName(29).name = 'R6';
LabelName(27).name = 'R7';
LabelName(3).name = 'R8';
LabelName(14).name = 'R9';
LabelName(38).name = 'R10';
LabelName(8).name = 'R11';
LabelName(30).name = 'R12';
LabelName(9).name = 'R13';


%% Plot the ranked, labeled mIBs
subplot(3,1,3)
for i_c = 1:Nc2
    delta = (i_c-1)*0.2;
    
    
    
    for i_b = 1:length(mIBs_Blue_Sorted)  % Put on Blue mIBs
        tempP = mIBs_Blue_Sorted(i_c+1,i_b);
        plot(delta + [i_b i_b], [0 -log10(tempP)], [clr(i_c) '-'],'linewidth',3)
        hold on
        if i_c == 3
            text(i_b,-1, sprintf('%2.0f',mIBs_Blue_Sorted(1,i_b)),'fontsize',6)
            mIB_num = mIBs_Blue_Sorted(1,i_b);
            try
                text(i_b, 1.5 - log10(min(mIBs_Blue_Sorted(2:4,i_b))), sprintf(LabelName(mIB_num).name) ,'fontsize',6,'color',[0 0 1])
            catch
            end
            text(i_b, 1 - log10(min(mIBs_Blue_Sorted(2:4,i_b))), 'p=','fontsize',6)
            text(i_b, 0.5 - log10(min(mIBs_Blue_Sorted(2:4,i_b))), sprintf('%1.0e', min(mIBs_Blue_Sorted(2:4,i_b)')),'fontsize',4)
        end
    end
    
    
    try
        [aG bG] = size(mIBs_Green_Sorted);
    catch
        bG = 0;  % there were no green channel most significant metrics
    end
    
    for i_b = 1:bG  % Put on Green mIBs
        tempP = mIBs_Green_Sorted(i_c+1,i_b);
        plot(length(mIBs_Blue_Sorted) + delta + [i_b i_b], [0 -log10(tempP)], [clr(i_c) '-'],'linewidth',3)
        hold on
        if i_c == 3
            text(length(mIBs_Blue_Sorted) + i_b,-1, sprintf('%2.0f',mIBs_Green_Sorted(1,i_b)),'fontsize',6)
            mIB_num = mIBs_Green_Sorted(1,i_b);
            try
                text(length(mIBs_Blue_Sorted) + i_b, 1.5 - log10(min(mIBs_Green_Sorted(2:4,i_b))), sprintf(LabelName(mIB_num).name) ,'fontsize',6,'color',[0 1 0])
            catch
            end
            text(length(mIBs_Blue_Sorted) + i_b, 1 - log10(min(mIBs_Green_Sorted(2:4,i_b))), 'p=','fontsize',6)
            text(length(mIBs_Blue_Sorted) + i_b, 0.5 - log10(min(mIBs_Green_Sorted(2:4,i_b))), sprintf('%1.0e', min(mIBs_Green_Sorted(2:4,i_b)')),'fontsize',4)
        end
    end
    
    try
        [aR bR] = size(mIBs_Red_Sorted);
    catch
        bR = 0;  % there were no red channel most significant metrics
    end
    
    for i_b = 1:bR % Put on Red mIBs
        i_b2 = length(mIBs_Red_Sorted) - i_b + 1;
        tempP = mIBs_Red_Sorted(i_c+1,i_b2);
        plot(length(mIBs_Blue_Sorted) + bG +delta + [i_b2 i_b2], [0 -log10(tempP)], [clr(i_c) '-'],'linewidth',3)
        hold on
        if i_c == 3
            text(length(mIBs_Blue_Sorted) + bG +i_b2, -1, sprintf('%2.0f',mIBs_Red_Sorted(1,i_b2)),'fontsize',6)
            mIB_num = mIBs_Red_Sorted(1,i_b2);
            try
                text(length(mIBs_Blue_Sorted) + bG + i_b2, 1.5 - log10(min(mIBs_Red_Sorted(2:4,i_b2))), sprintf(LabelName(mIB_num).name) ,'fontsize',6,'color',[1 0 0])
            catch
            end
            text(length(mIBs_Blue_Sorted) +bG +  i_b2, 1 - log10(min(mIBs_Red_Sorted(2:4,i_b2))), 'p=','fontsize',6)
            text(length(mIBs_Blue_Sorted) + bG + i_b2, 0.5 - log10(min(mIBs_Red_Sorted(2:4,i_b2))), sprintf('%1.0e', min(mIBs_Red_Sorted(2:4,i_b2)')),'fontsize',4)
        end
    end
    
    for i_b = 1:length(MCp_Sorted) % Put on Multicolor mIBs
        i_b2 = length(MCp_Sorted) - i_b + 1;
        tempP = MCp_Sorted(2,i_b2);
        plot(length(mIBs_Blue_Sorted) + bG + length(mIBs_Red_Sorted) +  [i_b2 i_b2] , [0 -log10(tempP)], [clr(i_c+1) '-'],'linewidth',6)
        mIB_num = MCp_Sorted(1,i_b2);
        try
            text(length(mIBs_Blue_Sorted) + bG + length(mIBs_Red_Sorted) + i_b2 -0.2, 1.5 - log10(tempP), sprintf(LabelName(mIB_num).name) ,'fontsize',6,'color',[1 0 1])
        catch
        end
        text(length(mIBs_Blue_Sorted) + bG + length(mIBs_Red_Sorted) +  i_b2 -0.2, 1 - log10(tempP), 'p=','fontsize',6)
        text(length(mIBs_Blue_Sorted) +bG + length(mIBs_Red_Sorted) +  i_b2 -0.2, 0.5 - log10(tempP), sprintf('%1.0e', tempP),'fontsize',6)
        text(length(mIBs_Blue_Sorted) + bG +length(mIBs_Red_Sorted) +  i_b2, -1, sprintf('%2.0f',MCp_Sorted(1,i_b2)),'fontsize',6)
    end
end

axis([0 (sum(good) + i_count_mc +1) 0 50])
plot([0 pos],[-log10(0.05) -log10(0.05)],'k--')
ylabel('-log10(p)')
xlabel('mIB Number, Ordered')

Max_SCmIBs = max(max(max(-log10(pval_out(pval_out>0)))));
Max_MCmIBs = max(max(  -log10(MCp(MCp>0)) ));

Grand_Max = max([Max_SCmIBs Max_MCmIBs]);

subplot(3,1,1)
axis([ 0 Nb2*1.1 0 Grand_Max*1.3  ])
subplot(3,1,2)
axis([ 0 (length(mIBs_Blue_Sorted) + bG +length(mIBs_Red_Sorted) +  i_b)*1.1 0 Grand_Max*1.3  ])
subplot(3,1,3)
axis([ 0 (length(mIBs_Blue_Sorted) + bG +length(mIBs_Red_Sorted) +  i_b) 0 Grand_Max*1.3  ])


disp([sprintf('%4.0f Lesions failed: ', length(n_nothing(n_nothing>0)) +  length(n_NotAll(n_NotAll>0)) + length(n_NaN(n_NaN>0)) )])
disp([sprintf('%4.0f Lesions with no mIB vals: ', length(n_nothing(n_nothing>0))) num2str(n_nothing)])
disp([sprintf('%4.0f Lesions with not all mIB vals: ', length(n_NotAll(n_NotAll>0))) num2str(n_NotAll)])
disp([sprintf('%4.0f Lesions with NaN mIB vals: ', length(n_NaN(n_NaN>0))) num2str(n_NaN)])

mIBs_Red_Sorted_stripped = mIBs_Red_Sorted;
mIBs_Green_Sorted_stripped = mIBs_Green_Sorted;
mIBs_Blue_Sorted_stripped = mIBs_Blue_Sorted;

redmin = min(mIBs_Red_Sorted, [], 1);
good = redmin <= .05;
mIBs_Red_Sorted_stripped = mIBs_Red_Sorted(:,good);
rs = size(mIBs_Red_Sorted_stripped(1,:),2);

greenmin = min(mIBs_Green_Sorted, [], 1);
good = greenmin <= .05;
mIBs_Green_Sorted_stripped = mIBs_Green_Sorted(:,good);
gs = size(mIBs_Green_Sorted_stripped(1,:),2);

bluemin = min(mIBs_Blue_Sorted, [], 1);
good = bluemin <= .05;
mIBs_Blue_Sorted_stripped = mIBs_Blue_Sorted(:,good);
bs = size(mIBs_Blue_Sorted_stripped(1,:),2);

MCmin = min(MCp_Sorted, [], 1);
good = MCmin <= .05;
mIBs_MC_Sorted_stripped = MCp_Sorted(:,good);
mcs = size(mIBs_MC_Sorted_stripped(1,:),2);

listmIBs(:,2) = [mIBs_Red_Sorted_stripped(1,:), mIBs_Green_Sorted_stripped(1,:), mIBs_Blue_Sorted_stripped(1,:),mIBs_MC_Sorted_stripped(1,:)]';
listmIBs(1:rs,1) = 1;
listmIBs(rs+1:rs+gs,1) = 2;
listmIBs(rs+gs+1:rs+gs+bs+mcs,1) = 3;

Matrix_Out3 = cat(3, Matrix_Out2(1,:,mIBs_Red_Sorted_stripped(1,:)), Matrix_Out2(2,:,mIBs_Green_Sorted_stripped(1,:)), ...
    Matrix_Out2(3,:,horzcat(mIBs_Blue_Sorted_stripped(1,:), mIBs_MC_Sorted_stripped(1,:))));

Matrix_Out3 = squeeze(Matrix_Out3);
list_mIBs = listmIBs';

cd Data
cd(folderDIR)

save Matrix_mIBs_3 Matrix_Out3 imID listmIBs LookType1 LookType2  n_nothing n_NotAll n_NaN
 OutFileName = ['Matrix_mIBs_3_' SaveName];
% statement = ['save ' OutFileName ' Matrix_Out3 imID listmIBs LookType1 LookType2 OutFileName n_nothing n_NotAll n_NaN'];
% eval(statement)
% statement = ['save ' OutFileName '_MIBsigID list_mIBs LookType1 LookType2 OutFileName'];
% eval(statement)

% each row of listmiB (color channel, mIB column number in Matrix_Out2) corresponds to a mIB column in Matrix_Out3
cd Track_mIBs
statement = ['save ' date '_' OutFileName ' Matrix_Out3 imID listmIBs LookType1 LookType2  n_nothing n_NotAll n_NaN'];
eval(statement)
cd ../../..

nmel = 0;
nnev = 0;
imID(imID==0) = [];
for i = 1:length(imID)
    if Diagnoses(imID(i),2) == 1
        nmel = nmel + 1;
    end
    
    if Diagnoses(imID(i),2) == 0
        nnev= nnev + 1;
    end
end


% N_mIBs_didnt_work = n_nothing + n_NotAll + n_NaN;
% N_autocrop_didnt_work = sum(ManualCheck);
subplot(3,1,1)
title([ folderDIR sprintf([ ', ' num2str(length(dirdir_Dat)) ' images-- ' ])  sprintf('  using: %4.0f Melanomas and %4.0f nevi w/ %4.0f excluded for bad mIBS',  nmel,nnev, N_mIBs_didnt_work ) '   ' date ],'Interpreter','none')

cd Data
cd(folderDIR)
if ~isdir('Track_Sig')
    mkdir('Track_Sig')
end
cd Track_Sig
saveas(gcf,['pVals_' folderDIR '_' date '.jpg'], 'jpg')
cd ..

[a b] = size(Matrix_Out3);
Matrix_2XL = zeros(a,b+1);
Matrix_2XL(:,2:b+1) = Matrix_Out3;

for i = 1:a
    Matrix_2XL(i,1) = Diagnoses(imID(i),2);
end

xlswrite('OutForJoel.xls',Matrix_2XL)
cd(homedir)
