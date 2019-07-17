% autocrop.m
% won't work for last 500 files of S2_MSKCC, runs out of memory
clear all
close all
homedir = pwd;


folderDIR = 'S0_Orig_Marghoob';
% started after clear border on Puig 614
range = [2]; % 3: S0_Orig_Marghoob [120]       2: S1_MSKCC   [900]      1: S2_Puig  [1081]  1: S3_DebbieMiller   5:S4b_NormNevi_RU_Mmx  1:S6_Pellacani    39:S11_UCI_OHS
macscreen =0 ; % SET ONLY ONE OF THE FOLLOWING VARIABLES TO "1"
PCscreen =0;
PCscreen2 = 0;
PCscreen_remote =0;

Smoothit = 0; % if curve fitting toolbox present
plotit = 1;
targetDIR = 'Data'; 
Hold_Off_Recrop = 0; % do not change
Tf = 1; % Threshold factor, for use in Run2
carryover = 0;
useblue = 0;

cd(targetDIR)              % Change Directory to folder where data lies
cd(folderDIR)
DataFolder = pwd;
load('GSdiag.mat')

Nims_bmp = length(dir('*.bmp'));
Nims_jpg = length(dir('*.jpg'));
Nims_tif = length(dir('*.tif'));

if Nims_bmp > Nims_jpg + Nims_tif
    Dat_Type = 'bmp';
elseif Nims_jpg > Nims_bmp + Nims_tif
    Dat_Type = 'jpg';
else
    Dat_Type = 'tiff';
end

dirdir_Dat = dir(sprintf(['*.' Dat_Type]));
cd(homedir)
C1 = 0.5;
C2 = 0.5; % constant for number of SD's +/- to make random noise in frame

save CurrentWorkingDir DataFolder dirdir_Dat C1 C2 folderDIR Dat_Type PCscreen PCscreen2 PCscreen_remote macscreen Tf carryover Smoothit useblue
cd(DataFolder)
Nims = length(dirdir_Dat);

flag_Nims =1;
first_time = 1;
num_done =0;

tic % needed for timed displayc

for i_im = range
    
    num_done = num_done + 1;
    
    try % mark if crop fails
        
        cd(DataFolder)
        name = dirdir_Dat(i_im).name;
        dat = imread(name);
        dat = double(dat);
        [Ny, Nx, Nz] = size(dat);
        cd(homedir)
        F0_fix_marghoob
        cd(homedir)
        F1_Init_Proc_Img
        
        F2_Re_Segment2;
        F7_write_ims
        Diagnoses(i_im,3) = 0;
        
        if first_time
            time_catch = toc;
            disp(sprintf('%5.0f[s] Now Autocropping image %3.0f, which is #%3.0f of %3.0f', time_catch,   i_im, num_done, (max(range) - min(range))  ))
        else
            time_catch = toc;
            exp_s = time_catch/num_done * ((max(range) - min(range) - num_done) + 1);
            exp_m = floor(exp_s/60);
            ss = exp_s-60*exp_m;
            exp_h = floor(exp_m/60);
            mm = floor(exp_s/60) - exp_h*60;
            disp(sprintf('%5.0f[s] Elapsed ... Now Autocropping image %3.0f, which is #%3.0f of %3.0f ... done in %2.0f:%2.0f:%2.0f [h:m:s]', time_catch,   i_im, num_done, (max(range) - min(range)), exp_h,mm,ss))
        end
        Diagnoses(i_im,9) = 0; % Save that there was a mIB code success
    catch
        
        Diagnoses(i_im,9) = 1; % Save that there was a mIB code failure
        Diagnoses(i_im,20) = 1; % Save that there was a mIB code failure
        disp(sprintf('AUTO-CROP FAIL ON IMAGE %5.0f',i_im))
        
   end
   
    
    first_time = 0;
    cd(homedir)
    cd Data
    cd(folderDIR)
    save GSdiag Diagnoses   
    cd ../..
end

cd(homedir)

close all
beep; pause(0.2);beep; pause(0.2);beep; pause(0.2);
disp(['Done Auto-Cropping data set: ' folderDIR])


