%%Clean_Up%%%   Wiggling Fingers at the Force-Length Relationship: 
clearvars   %   FDI, Architectural Influence on Force Potential.
clc         %   Zackary Scalyer
close all   %   The Pennsylvania State University - Berks
%%%%%%%%%%%%%   Spring_17: Dr. Ben Infantolino
%{
    Purpose - To filter EMG and Force data of the FDI, and find ultrasound
        frames corresponding to 20 percent increments of Force and EMG. 
        
    input
    -------
    ultrasound video of muscle fascicles (avi)
    EMG & force voltage data (csv)

    output
    --------
    Processing - 3x3 Figure of data processing stages for EMG and force
             - Five Graphs for each EMG signal processing 
               1) Raw data (EMG)                                
               2) High-pass data (20hz cut-off)                 
               3) Root Mean Square in 0.1s window                        
               4) Smoothed (1 Hz cut-off)                       
               5) Normalized (with respect to max of each Ramp)
             - Four graphs for each Force signal processing   
               6) Raw data                                     
               7) Low-pass data (1hz cut-off)                   
               8) Magnitude of LP_Force (subtrat min-value)     
               9) Normalized (with respect to max of each Ramp)
    FData - Data set of filter data for each joint angle
        - Column 1: Time (Time)
        - Column 2: Normalised EMG (NormEMG)
        - Column 3: Normalised Force (NormForce)
        - Column 4: Magnitude Force (MagForce)
    A_BAD - Data set of sub-max activation percentage for each joint angle
        - Column 1: Time
        - Column 2: Normalised EMG
        - Column 3: Normalised Force
        - Column 4: Magnitude Force
    F_BAD - Data set of sub-max force percentage for each joint angle
        - Column 1: Time
        - Column 2: Normalised EMG
        - Column 3: Normalised Force
        - Column 4: Magnitude Force
    SubMaxPlot - 1x2 Figure of force & EMG with sub-max percentages 
             - 1) text is percent force at percnet activation  
             - 2) text is percent activation at percnet force

    Notes
    -------
    i. filtmat_class.m is a Butterworth filter used as a subroutine and was 
     written by Dr. John Challis (1997)



    Updates
    i. 7/22/2017: Changed EMG filtering, RMS insed of abs.
    ---------
%}
%% Initialize variables
%-------------------------%
% Enter hard drive number %
%-------------------------%
Drive = 'E';
%----------------------%
% Enter subject number %
%----------------------%
n = [3];
%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st MCP joint angles %%
%%%%%%%%%%%%%%%%%%%%%%%%%
MCP_drg = [0 5 10 15 20];
%%%%%%%%%%%%%%%%%%%%%
% Ramp contraction %%
%%%%%%%%%%%%%%%%%%%%%
Ramp = [1 2];
%%%%%%%%%%%%%%%%%%%%
% Sub-max percent %%
%%%%%%%%%%%%%%%%%%%%
Percent = 0.2;

disp('Finished Initialize variables');
%% Add path
%---------------------------%
% Enter subject folder path %
%---------------------------%
for i = 1:length(n)
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d',Drive, n(i)));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d\\n%d_out',Drive, n(i), n(i)));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d\\n%d_filtered',Drive, n(i), n(i)));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\Calibration',Drive));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Library',Drive));
end
clear i


disp('Finished Add path');

%% Load Cell Calibration
load('WigglyLoad.mat');

disp('Finished loading load cell calibration');
%% Read Subject's data
%%%%%%%%%%%%%%%%%%%%%%%%%
% CSV of EMG and Force %%
%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        for j = 1:length(Ramp)
            fname = sprintf('S%d_ZS_%d_%d.csv',n(I),MCP_drg(i),Ramp(j));
            eval(sprintf('S%d_%d_%d = csvread(fname,40,0);',n(I),MCP_drg(i),Ramp(j)));
        end
    end
end
clear I i j fname;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample Rate for EMG and Force data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(sprintf('dt = S%d_20_1(2,1)- S%d_20_1(1,1);',n(1),n(1)));
eval(sprintf('time = S%d_20_1(:,1);',n(1),n(1)));
%%%%%%%%%%%%%%%%%%%
% avi ultrasound %%
%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        for j = 1:length(Ramp)
            fname_m = sprintf('S%d_ZS_%d_%d.avi',n(I),MCP_drg(i),Ramp(j));
            eval(sprintf('US_S%d_%d_%d = VideoReader(fname_m);',n(I),MCP_drg(i),Ramp(j))); % US Video
            eval(sprintf('num_US_S%d_%d_%d = US_S%d_%d_%d.NumberOfFrames;',...          % number of frames
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('FR_US_S%d_%d_%d = US_S%d_%d_%d.FrameRate;',...                % frame rate
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('dt_FR_US_S%d_%d_%d = 1/FR_US_S%d_%d_%d;',...                  % time interval
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('Time_US_S%d_%d_%d = 0:dt_FR_US_S%d_%d_%d:(num_US_S%d_%d_%d*dt_FR_US_S%d_%d_%d)-dt_FR_US_S%d_%d_%d;',... % time array for each video
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));       
        end
    end
end
clear I i j fname_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove module level variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        for j = 1:length(Ramp)
            eval(sprintf('clear FR_US_S%d_%d_%d;', n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('clear dt_FR_US_S%d_%d_%d;', n,MCP_drg(i),Ramp(j)));
        end
    end
end
clear I i j;

disp('Finished Read Subject''s data');
%% Truncate data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncate to length of US avi %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        for j = 1:length(Ramp)
            eval(sprintf('[x] = find(S%d_%d_%d(:,1) >= max(Time_US_S%d_%d_%d),1);',... % time US ended
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j))); % 
            x=x-1; % reduce time by one sample b/c last US frame my have already occured
            eval(sprintf('Data_S%d_%d_%d = zeros(x,4);',... % initialise empty data set
                n(I),MCP_drg(i),Ramp(j))); 
            eval(sprintf('Data_S%d_%d_%d(:,1) = S%d_%d_%d(1:x,1);',... % truncate time
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('Data_S%d_%d_%d(:,2) = S%d_%d_%d(1:x,2);',... % truncate EMG
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('Data_S%d_%d_%d(:,3) = S%d_%d_%d(1:x,10);',... % truncate Force
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('Data_S%d_%d_%d(:,4) = V2F(Data_S%d_%d_%d(:,3));',... % convert voltage to force (N)
                n(I),MCP_drg(i),Ramp(j),n(I),MCP_drg(i),Ramp(j)));
        end
    end
end
clear I i j x;
%-----------------------------------------------%
% Select ramp trial with highest absolute force %
%-----------------------------------------------%
for r = 1:length(n)
    for i = 1:length(MCP_drg)
        reply = 1;
        for j = 1:length(Ramp) % Modual level low pass force.
            eval(sprintf('[LP_Force_S%d_%d_%d] = filtmat_class(dt,1,Data_S%d_%d_%d(:,4),1);',...
                n(r),MCP_drg(i),Ramp(j),n(r),MCP_drg(i),Ramp(j)));
        end
        eval(sprintf('numOfFram1_m = num_US_S%d_%d_1;',n(r),MCP_drg(i))); 
        eval(sprintf('img1_m = US_S%d_%d_1;',n(r),MCP_drg(i)));
        eval(sprintf('numOfFram2_m = num_US_S%d_%d_2;',n(r),MCP_drg(i))); 
        eval(sprintf('img2_m = US_S%d_%d_2;',n(r),MCP_drg(i)));
        figure(1);
        ax(1) = subplot(2,2,1);
        for I = 1:numOfFram1_m
           img1_read_m = read(img1_m,I); 
           imshow(img1_read_m,[])
        end
        title(['Slect to replay AVI for: S',int2str(n(r)),' drg', int2str(MCP_drg(i)),...
        ' 1st Ramp']);
        ax(2) = subplot(2,2,2);
        for J = 1:numOfFram2_m
           img2_read_m = read(img2_m,J); 
           imshow(img2_read_m,[])
        end
        title(['Slect to replay AVI for: S',int2str(n(r)),' drg', int2str(MCP_drg(i)),...
        ' 2nd Ramp']);
        ax(3) = subplot(2,2,3);
        eval(sprintf('plot(Data_S%d_%d_1(:,1),LP_Force_S%d_%d_1);',n(r),MCP_drg(i),n(r),MCP_drg(i)));
        title(['1Hz LowPass Force: S',int2str(n(r)),' drg', int2str(MCP_drg(i)),...
        ' 1st Ramp']);
        ylabel('Force(N)');
        xlabel('Time(S)');
        ax(4) = subplot(2,2,4);
        eval(sprintf('plot(Data_S%d_%d_2(:,1),LP_Force_S%d_%d_2);',n(r),MCP_drg(i),n(r),MCP_drg(i)));
        title(['1Hz LowPass Force: S',int2str(n(r)),' drg', int2str(MCP_drg(i)),...
        ' 2nd Ramp']);
        ylabel('Force(N)');
        xlabel('Time(S)');
        linkaxes([ax(3),ax(4)],'y');
        while reply == 1
            [x] = ginput(1);
            axesSlect = gca;
            chooseData = (ax == axesSlect);
            if chooseData(1)== 1
                ax(1) = subplot(2,2,1);
                for ii = 1:numOfFram1_m
                    img1_read_m = read(img1_m,ii); 
                    imshow(img1_read_m,[])
                end
                title(['Slect to replay AVI for: S',int2str(n(r)),' drg', int2str(MCP_drg(i)),...
                    ' 1st Ramp']);
                reply = 1;
            end
            if chooseData(2)== 1
                ax(2) = subplot(2,2,2);
                for jj = 1:numOfFram2_m
                    img2_read_m = read(img2_m,jj); 
                    imshow(img2_read_m,[])
                end
                title(['Slect to replay AVI for: S',int2str(n(r)),' drg', int2str(MCP_drg(i)),...
                    ' 2nd Ramp']);
                reply = 1;
            end
            if chooseData(3)== 1
                eval(sprintf('Data_S%d_%d = Data_S%d_%d_1;',...
                    n(r),MCP_drg(i),n(r),MCP_drg(i)));
                eval(sprintf('US_S%d_%d = US_S%d_%d_1;',...
                    n(r),MCP_drg(i),n(r),MCP_drg(i)));
                eval(sprintf('Time_US_S%d_%d = Time_US_S%d_%d_1;',...
                    n(r),MCP_drg(i),n(r),MCP_drg(i)));
                reply = 0;
                % save a copy of data trial that will be analised
                path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                    int2str(n(r)),'\S',int2str(n(r)),'_ZS_',int2str(MCP_drg(i)),'_ramp_1.csv'];
                dat = ['Data_S',int2str(n(r)),'_',int2str(MCP_drg(i))];
                csvwrite(path,dat); 
            end
            if chooseData(4)== 1
                eval(sprintf('Data_S%d_%d = Data_S%d_%d_2;',...
                    n(r),MCP_drg(i),n(r),MCP_drg(i)));
                eval(sprintf('US_S%d_%d = US_S%d_%d_2;',...
                    n(r),MCP_drg(i),n(r),MCP_drg(i)));
                eval(sprintf('Time_US_S%d_%d = Time_US_S%d_%d_2;',...
                    n(r),MCP_drg(i),n(r),MCP_drg(i)));
                reply = 0;
                % save a copy of data trial that will be analised
                path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                    int2str(n(r)),'\S',int2str(n(r)),'_ZS_',int2str(MCP_drg(i)),'_ramp_2.csv'];
                dat = ['Data_S',int2str(n(r)),'_',int2str(MCP_drg(i))];
                csvwrite(path,dat); 
            end
        end
        close 1;
    end
end
clear r i j I J ii jj reply x ax clickedax chooseData...
    img1_m img2_m img1_read_m img2_read_m path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove module level variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        for j = 1:length(Ramp)
            eval(sprintf('clear LP_Force_S%d_%d_%d;', n(I),MCP_drg(i),Ramp(j)));
            eval(sprintf('clear num_US_S%d_%d_%d;', n,MCP_drg(i),Ramp(j)));
        end
    end
end
clear I i j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove raw data from enviorment %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        for j = 1:length(Ramp)
            eval(sprintf('clear S%d_%d_%d;', n(I),MCP_drg(i),Ramp(j)));
        end
    end
end
clear I i j;

disp('Finished Truncate data');
%% Filter EMG Data
%   i) High-pass data (20Hz cut-off)
%   ii) Root Mean Square in 0.1s window  
%   iii) Smoothed (1 Hz LowPass cut-off)
%   iv) Normalized (with respect to maximum of each Ramp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20Hz High-pass filter raw EMG data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('[HP_EMG_S%d_%d] = filtmat_class(dt,20,Data_S%d_%d(:,2),2);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Root Mean Square in 0.1s window %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        dat = eval(sprintf('HP_EMG_S%d_%d',n(I),MCP_drg(i)));    % grab EMG data
        tim = eval(sprintf('Data_S%d_%d(:,1)',n(I),MCP_drg(i))); % grab EMG time
        windowEnd = find(tim >= 0.1, 1);                   % window of 0.1s
        append_dat = [zeros(floor(.5*windowEnd),1);...     % insirt zeros of 1/2 indow length to the frount of the EMG singla
            dat;...                                        % insert EMG signal
            zeros(floor(.5*windowEnd),1)];                 % add 1/2 of the window length to the end of EMG singnal
        rms_dat = zeros(length(append_dat),1);             % initsalise vetor of zeros    
        for j = 1:length(rms_dat)                          % windowed RMS
            if j+windowEnd <= length(rms_dat)
                rms_dat(j) = rms(append_dat(j:j+windowEnd,1));
            end
        end
        timShift = tim(floor(.5*windowEnd):end);          % adjusted time 
        rms_dat = rms_dat(1:length(timShift));            % truncate EMG to length of time
        eval(sprintf('RMS_HP_EMG_S%d_%d = rms_dat;',...
            n(I),MCP_drg(i)));                            % Name subject joint angle EMG signal
        eval(sprintf('Time_S%d_%d = Data_S%d_%d(:,1);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));            % Un-adjusted time for plots
         eval(sprintf('Data_S%d_%d = Data_S%d_%d(floor(.5*windowEnd):end,:);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));            % Truncate original data to EMG signal
    end
end
clear I i j dat tim windowEnd append_dat rms_dat timShift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoothing/1Hz LowPass of R_H__EMG data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('[LP_RMS_HP_EMG_S%d_%d] = filtmat_class(dt,1,RMS_HP_EMG_S%d_%d,1);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of EMG data         %%
% with respect to max of each Ramp  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('Norm_LP_RMS_HP_EMG_S%d_%d = LP_RMS_HP_EMG_S%d_%d/max(LP_RMS_HP_EMG_S%d_%d);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;

disp('Finished Filter EMG Data');
%% Filter Force Data
%   i) 1Hz low-pass filter force
%   ii) Magnitude of LP_Force (Subtration of min-value)
%   iii) Normalized (with respect to maximum of each Ramp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1Hz low-pass filter force data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('[LP_Force_S%d_%d] = filtmat_class(dt,1,Data_S%d_%d(:,3),1);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnitude of LP_Force data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('Min_LP_Force_S%d_%d = min(LP_Force_S%d_%d);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
        eval(sprintf('M_LP_Force_S%d_%d = LP_Force_S%d_%d - Min_LP_Force_S%d_%d;',...
            n(I),MCP_drg(i),n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of Force data      %%
% with respect to max of each Ramp %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('Max_M_LP_Force_S%d_%d = max(M_LP_Force_S%d_%d);',...
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
        eval(sprintf('Norm_M_LP_Force_S%d_%d = M_LP_Force_S%d_%d/Max_M_LP_Force_S%d_%d;',...
            n(I),MCP_drg(i),n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;

disp('Finished Filter Force Data');
%% 3x3 Figrure of Data Processing 
% Produece 5 Graphs for each EMG signal processing 
%  1) Raw data (EMG)                                
%  2) High-pass data (20hz cut-off)                 
%  3) Root Mean Square in 0.1s window                     
%  4) Smoothed (1 Hz cut-off)                       
%  5) Normalized (with respect to max of each Ramp)
% With 4 graphs for each Force signal processing   
%  6) Raw data                                     
%  7) Low-pass data (1hz cut-off)                   
%  8) Magnitude of LP_Force (subtrat min-value)     
%  9) Normalized (with respect to max of each Ramp) 
k=6; 
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            l = 1;
            figure(k);
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), Data_S%d_%d(:,2));',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['Raw EMG: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('EMG(Voltge)');
            grid on; 
            l = l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Time_S%d_%d, HP_EMG_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['20Hz HP EMG: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('EMG(Voltge)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), RMS_HP_EMG_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['0.1s windowed RMS: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('EMG(Voltge)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), LP_RMS_HP_EMG_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['Smoothed 1Hz LP EMG: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('EMG(Voltge)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), Norm_LP_RMS_HP_EMG_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['Normalize EMG: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('EMG(Voltge)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), Data_S%d_%d(:,4));',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['Raw Force: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('Force(N)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), LP_Force_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['1Hz LP Force: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('Force(N)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), M_LP_Force_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['Magnitude of LP Force: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('Force(N)');
            grid on; 
            l=l+1;
            subplot(3,3,l);
            eval(sprintf('plot(Data_S%d_%d(:,1), Norm_M_LP_Force_S%d_%d);',...
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            title(['Normalized Force: ', int2str(MCP_drg(i)),' drg']);
            xlabel('Time(S)');
            ylabel('Force(N)');
            grid on; 
            k = k+1;

    end
end
clear I k l i;
%%%%%%%%%%%%%%%%%%%%%%
% Save Figures 6:15 %%
%%%%%%%%%%%%%%%%%%%%%%
k = 6;
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                int2str(n(I)),'\n',int2str(n(I)),'_filtered\Processing_n',...
                int2str(n(I)),'_',int2str(MCP_drg(i)),'drg.fig'];
            eval(sprintf('saveas(%d,path);',k));
            k = k+1;
    end
end
clear I i k path;

disp('Finished 3x3 Figrure of Data Processing');
%% Create set of filtered data
%{
  FData - Data set of filter data for each Ramp at each joint angle
        - Column 1: Time
        - Column 2: Raw EMG
        - Column 3: Raw Force
        - Column 4: Normalised EMG
        - Column 6: Normalised Force
        - Column 7: Absolute Force
%}
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('FData_S%d_%d = zeros(size(Data_S%d_%d));',... % initialise empty data set
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
        eval(sprintf('FData_S%d_%d(:,1) = Data_S%d_%d(:,1);',... % Time
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
        eval(sprintf('FData_S%d_%d(:,2) = Norm_LP_RMS_HP_EMG_S%d_%d;',... % Normalised EMG
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
        eval(sprintf('FData_S%d_%d(:,3) = Norm_M_LP_Force_S%d_%d;',... % Normalised Force
            n(I),MCP_drg(i),n(I),MCP_drg(i)));
        eval(sprintf('FData_S%d_%d(:,4) = M_LP_Force_S%d_%d;',... % Magnitude Force
        n(I),MCP_drg(i),n(I),MCP_drg(i)));
    end
end
clear I i;
%%%%%%%%%%%%%%%%%%%%%%%
% Save filtered data %%
%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                int2str(n(I)),'\n',int2str(n(I)),'_filtered\FData_n',int2str(n(I)),'_',...
                int2str(MCP_drg(i)),'drg.csv'];
            eval(sprintf('Time = Data_S%d_%d(:,1);',n(I),MCP_drg(i))); % Time
            eval(sprintf('RawEMG = Data_S%d_%d(:,2);',n(I),MCP_drg(i))); % Raw EMG
            eval(sprintf('RawForce =Data_S%d_%d(:,3);',n(I),MCP_drg(i))); % Raw Force
            eval(sprintf('NormEMG = Norm_LP_RMS_HP_EMG_S%d_%d;',n(I),MCP_drg(i))); % Normal EMG
            eval(sprintf('NormForce = Norm_M_LP_Force_S%d_%d;',n(I),MCP_drg(i))); % Normal Force
            eval(sprintf('MagForce = M_LP_Force_S%d_%d;',n(I),MCP_drg(i))); % Magnitude Force
            t = table(Time,RawEMG,RawForce,NormEMG,NormForce,MagForce);
            writetable(t,path);
    end
end
clear I i path t Time NormEMG NormForce MagForce;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove raw and un-used filted data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('clear Data_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear HP_EMG_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear R_HP_EMG_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear LP_R_HP_EMG_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear Max_LP_R_HP_EMG_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear Norm_LP_R_HP_EMG_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear LP_Force_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear M_LP_Force_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear Min_LP_Force_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear Max_M_LP_Force_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear Norm_M_LP_Force_S%d_%d;', n(I),MCP_drg(i)));
    end
end
clear I i;

disp('Finished Create set of filtered data');
%% Find Sub-max force by percent activation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find sub-max levels of activation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        clear EmgVector;
        k = 0;          % counter
        hitMaxEMG = 0;  % logical
        ii = 2;         % index to compare to
        EmgVector = eval(sprintf('FData_S%d_%d(:,2)',...
            n(I),MCP_drg(i))); % grab emg data
        while hitMaxEMG == 0 
            if EmgVector(ii) == max(EmgVector)  % when EMG hits maximum
                hitMaxEMG = 1;                  % stop looking for percentages of EMG
                k = k+1;    % count number of times EMG incred by 20 percent
                eval(sprintf('A_BAD_%d(%d,1:4) = FData_S%d_%d(%d,1:4);',... % grap FData 
                    MCP_drg(i),k,n(I),MCP_drg(i),ii));
            elseif floor(EmgVector(ii)/Percent) > floor(EmgVector(ii-1)/Percent) % when EMG increses by 20 percent
                k = k+1;    % count number of times EMG incred by percent
                eval(sprintf('A_BAD_%d(%d,1:4) = FData_S%d_%d(%d,1:4);',... % grap FData 
                    MCP_drg(i),k,n(I),MCP_drg(i),ii));
            end
            ii = ii+1;  % compare to next index
        end
    end
end

clear I i k hitMaxForce EmgVector

%%%%%%%%%%%%%%%%%%%%%%
% Check SubAct Data %%
%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(MCP_drg)
    dat = eval(sprintf('A_BAD_%d',MCP_drg(i))); % grab sub activation data
    if length(dat) ~= 5 % if nuber of EMG percentages found is not 5
        dat_inter = floor(dat(:,2).*10);     % sub activation as integers
        expectPercent = floor((Percent:Percent:1).*10);   % expected result
        for j = 1:length(expectPercent)
            count = sum(dat_inter == expectPercent(j)); % count frequencys
            if count > 1 % if sub activation occurs more than one time
                toMany = find(dat(:,2) >= (expectPercent(j)./10) &...
                    dat(:,2) < (expectPercent(j)./10+0.1)); % how may times?
                medTim = median(dat(toMany,1));                   % find the median time of all occures
                keepFirst = sort(abs(dat(toMany,1)-medTim));      % sort distence from median 
                remove = find(round(dat(:,1),4) == round(keepFirst(2:end)+medTim,4)); % find the time of all but first
                dat(remove,:) = [];                               % remove all but first
                eval(sprintf('A_BAD_%d = dat;',MCP_drg(i))); % re-save ended sub EMG data
            end
        end
    end
end
% clear i dat dat_inter expectPercent count toMany medTim keepFirst remove

%%%%%%%%%%%%%%%%%%%%%
% Save SubAct data %%
%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                int2str(n(I)),'\n',int2str(n(I)),'_filtered\SubAct_n',int2str(n(I)),'_',...
                int2str(MCP_drg(i)),'drg.csv'];
            eval(sprintf('Time = A_BAD_%d(:,1);',MCP_drg(i)));
            eval(sprintf('SubAct = A_BAD_%d(:,2);',MCP_drg(i)));
            eval(sprintf('SubForce = A_BAD_%d(:,3);',MCP_drg(i)));
            eval(sprintf('MagForce = A_BAD_%d(:,4);',MCP_drg(i)));
            t = table(Time,SubAct,SubForce,MagForce);
            writetable(t,path);
    end
end
clear I i path t Time SubAct SubForce MagForce;

disp('Finished Find Sub-max force by percent activation');
%% Find sub-max force by percent force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find sub-max levels of force %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        clear forceVector;
        k = 0;            % counter
        hitMaxForce = 0;  % logical
        ii = 2;           % index to compare to
        forceVector = eval(sprintf('FData_S%d_%d(:,3)',...
            n(I),MCP_drg(i))); % grab force data
        while hitMaxForce == 0
            if forceVector(ii) == max(forceVector) % when force hits maximum
                hitMaxForce = 1;                   % stop looking for percentages of force
                k = k+1;    % count number of times EMG incred by 20 percent
                eval(sprintf('F_BAD_%d(%d,1:4) = FData_S%d_%d(%d,1:4);',... % grap FData 
                    MCP_drg(i),k,n(I),MCP_drg(i),ii));
            elseif floor(forceVector(ii)/Percent) > floor(forceVector(ii-1)/Percent) % when force increses by 20 percent
                k = k+1;    % count number of times force incred by percent
                eval(sprintf('F_BAD_%d(%d,1:4) = FData_S%d_%d(%d,1:4);',...
                    MCP_drg(i),k,n(I),MCP_drg(i),ii));
            end 
            ii = ii+1;  % compare to next index
       end
    end
end
clear I i k hitMaxForce forceVector;
%%%%%%%%%%%%%%%%%%%%%%
% Check SubAct Data %%
%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(MCP_drg)
    dat = eval(sprintf('F_BAD_%d',MCP_drg(i))); % grab sub force data
    if length(dat) ~= 5 % if nuber of force percentages found is not 5
        dat_inter = floor(dat(:,2).*10);     % sub force as integers
        expectPercent = floor((Percent:Percent:1).*10);   % expected result
        for j = 1:length(expectPercent)
            count = sum(dat_inter == expectPercent(j)); % count frequencys
            if count > 1 % if force activation occurs more than one time
                toMany = find(dat(:,2) >= (expectPercent(j)./10) &...
                    dat(:,2) < (expectPercent(j)./10+0.1)); % how may times?
                medTim = median(dat(toMany,1));                   % find the median time of all occures
                keepFirst = sort(abs(dat(toMany,1)-medTim));      % sort distence from median 
                remove = find(round(dat(:,1),4) == round(keepFirst(2:end)+medTim,4)); % find the time of all but first
                dat(remove,:) = [];                               % remove all but first
                eval(sprintf('F_BAD_%d = dat;',MCP_drg(i))); % re-save ended sub force data
            end
        end
    end
end
% clear i dat dat_inter expectPercent count toMany medTim keepFirst remove
%%%%%%%%%%%%%%%%%%%%%%%
% Save SubForce data %%
%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                int2str(n(I)),'\n',int2str(n(I)),'_filtered\SubForce_n',int2str(n(I)),'_',...
                int2str(MCP_drg(i)),'drg.csv'];
            eval(sprintf('Time = F_BAD_%d(:,1);',MCP_drg(i)));
            eval(sprintf('SubAct = F_BAD_%d(:,2);',MCP_drg(i)));
            eval(sprintf('SubForce = F_BAD_%d(:,3);',MCP_drg(i)));
            eval(sprintf('MagForce = F_BAD_%d(:,4);',MCP_drg(i)));
            t = table(Time,SubAct,SubForce,MagForce);
            writetable(t,path);
    end
end
clear I i path t Time SubAct SubForce MagForce;

disp('Finished Find sub-max force by percent force');
