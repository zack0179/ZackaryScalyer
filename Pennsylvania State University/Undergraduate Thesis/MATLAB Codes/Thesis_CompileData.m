%%Clean_Up%%%   Wiggling Fingers at the Force-Length Relationship: 
clearvars   %   FDI, Architectural Influence on Force Potential.
clc         %   Zackary Scalyer
close all   %   The Pennsylvania State University - Berks
%%%%%%%%%%%%%   Summer_17: Dr. Ben Infantolino
%{
    Purpose - To compile and evaluatesubject arcatrue data inorder to make
        a forece length grap at percetanges of forec and activation 

    input
    -------
    BAD - sujects data set of all sub-max percentage for each joint angle
      - Column 1: Joint angle (Drg)
      - Column 2: Time (Time)
      - Column 3: Normalised EMG (NormEMG)
      - Column 4: Normalised Force (NormForce)
      - Column 5: Magnitude Force (MagForce)
      - Column 6: Muscle fascicle length (FacLength)
      - Column 7: Pennation angle (FacAngle)
      - Column 8: Muscle thickness (MusThick)
    

    output
    --------
    One figure for each subject evaluating sub-max perameters over fascicle
        length.
    A forec-length graph at sub-max force and activation perecntages for
        each subject and one for the average of all subjects. 

    One struct of all data

  Notes
  -------
    Open MATLAB current folder to the 'CompileData' folder containing this
        script and analysed subject data.
    V2 - 8-1-17; grabs non normlised data 
    V3 - 8-1-17; graphs are now normilised 
    V4 - 8-1-17; figures not saving, removed. LFopt and norm Lfopt now saved in struct for
        activation and force. 
%}
%% Establish Perameters

%%%%%%%%%%%%%%%%%
% Investigator %%
%%%%%%%%%%%%%%%%%
k = 5;
%%%%%%%%%%%%%
% Subjects %%
%%%%%%%%%%%%%
n = [1 3 4 5 6 8 11 12 16];


disp('Finished Cell: Establish Perameters');
%% Load Cell Calibration
load('WigglyLoad.mat');

disp('Finished loading load cell calibration');
%% A struct
Wiggles_4 = struct('Subject',[], 'Investigator', [], 'subAct_0', [], 'facLen_act_0', [], ... 
    'subFor_0', [],'facLen_for_0', [], 'subAct_5', [], 'facLen_act_5', [], ... 
    'subFor_5', [],'facLen_for_5', [],'subAct_10', [], 'facLen_act_10', [], ... 
    'subFor_10', [],'facLen_for_10', [],'subAct_15', [], 'facLen_act_15', [], ... 
    'subFor_15', [],'facLen_for_15', [],'subAct_20', [], 'facLen_act_20', [], ... 
    'subFor_20', [],'facLen_for_20', [], 'A_20', [], 'A_40',...
    [], 'A_60', [], 'A_80', [], 'A_100', [], 'Lf_A_20', [], 'Lf_A_40',...
    [], 'Lf_A_60', [], 'Lf_A_80', [], 'Lf_A_100', [], 'F_20', [], 'F_40',...
    [], 'F_60', [], 'F_80', [], 'F_100', [], 'Lf_F_20', [], 'Lf_F_40',...
    [], 'Lf_F_60', [], 'Lf_F_80', [], 'Lf_F_100', [], 'Lfopt_A_100', [],...
    'Lfopt_A_80', [], 'Lfopt_A_60', [], 'Lfopt_A_40', [], 'Lfopt_A_20', [],...
    'Lfopt_F_100', [],'Lfopt_F_80', [], 'Lfopt_F_60', [], 'Lfopt_F_40', [],...
    'Lfopt_F_20', [], 'norm_lf_A_100', [],'norm_lf_A_80', [],'norm_lf_A_60', [],...
    'norm_lf_A_40', [],'norm_lf_A_20', [],'norm_lf_F_100', [],'norm_lf_F_80', [],...
    'norm_lf_F_60', [],'norm_lf_F_40', [],'norm_lf_F_20', [], 'norm_FL_A_100', [],...
    'norm_FL_A_80', [],'norm_FL_A_60', [],'norm_FL_A_40', [],'norm_FL_A_20', [],...
    'norm_FL_F_100', [],'norm_FL_F_80', [],'norm_FL_F_60', [],'norm_FL_F_40', [],...
    'norm_FL_F_20', []);
%% Figures for subject data

for I = 1:length(n)
    % grab subject data
    file = ['BAD_n',num2str(n(I)),'k',num2str(k),'.csv']; % File name
    dat = csvread(file,1,0); % subject data
    
    Wiggles_4(I).Subject = n(I); % save subject nuber in struct
    Wiggles_4(I).Investigator = k;
    
    % get indice of each joint angle
    drg0 = find(dat(:,1)==0);       % MCP drg 0
    drg5 = find(dat(:,1)==5);   % MCP drg 5
    drg10 = find(dat(:,1)==10);  % MCP drg 10
    drg15 = find(dat(:,1)==15);  % MCP drg 15
    drg20 = find(dat(:,1)==20);  % MCP drg 20
    
    % sub-max data for each joint angle
    subAct_0 = dat(drg0(1:5),3);
    facLen_act_0 = dat(drg0(1:5),6);
    subFor_0 = dat(drg0(6:10),4);
    facLen_for_0 = dat(drg0(6:10),6);
    
    Wiggles_4(I).subAct_0 = dat(drg0(1:5),[2,3]); % save sub-max data in struct
    Wiggles_4(I).facLen_act_0 = dat(drg0(1:5),[2,6]);
    Wiggles_4(I).subFor_0 = dat(drg0(6:10),[2,4]);
    Wiggles_4(I).facLen_for_0 = dat(drg0(6:10),[2,6]);
    
    subAct_5 = dat(drg5(1:5),3);
    facLen_act_5 = dat(drg5(1:5),6);
    subFor_5 = dat(drg5(6:10),4);
    facLen_for_5 = dat(drg5(6:10),6);
    
    Wiggles_4(I).subAct_5 = dat(drg5(1:5),[2,3]); % save sub-max data in struct
    Wiggles_4(I).facLen_act_5 = dat(drg5(1:5),[2,6]);
    Wiggles_4(I).subFor_5 = dat(drg5(6:10),[2,4]);
    Wiggles_4(I).facLen_for_5 = dat(drg5(6:10),[2,6]);
    
    subAct_10 = dat(drg10(1:5),3);
    facLen_act_10 = dat(drg10(1:5),6);
    subFor_10 = dat(drg10(6:10),4);
    facLen_for_10 = dat(drg10(6:10),6);
    
    Wiggles_4(I).subAct_10 = dat(drg10(1:5),[2,3]); % save sub-max data in struct
    Wiggles_4(I).facLen_act_10 = dat(drg10(1:5),[2,6]);
    Wiggles_4(I).subFor_10 = dat(drg10(6:10),[2,4]);
    Wiggles_4(I).facLen_for_10 = dat(drg10(6:10),[2,6]);
    
    subAct_15 = dat(drg15(1:5),3);
    facLen_act_15 = dat(drg15(1:5),6);
    subFor_15 = dat(drg15(6:10),4);
    facLen_for_15 = dat(drg15(6:10),6);
    
    Wiggles_4(I).subAct_15 = dat(drg15(1:5),[2,3]); % save sub-max data in struct
    Wiggles_4(I).facLen_act_15 = dat(drg15(1:5),[2,6]);
    Wiggles_4(I).subFor_15 = dat(drg15(6:10),[2,4]);
    Wiggles_4(I).facLen_for_15 = dat(drg15(6:10),[2,6]);
    
    subAct_20 = dat(drg20(1:5),3);
    facLen_act_20 = dat(drg20(1:5),6);
    subFor_20 = dat(drg20(6:10),4);
    facLen_for_20 = dat(drg20(6:10),6);
  
    Wiggles_4(I).subAct_20 = dat(drg20(1:5),[2,3]); % save sub-max data in struct
    Wiggles_4(I).facLen_act_20 = dat(drg20(1:5),[2,6]);
    Wiggles_4(I).subFor_20 = dat(drg20(6:10),[2,4]);
    Wiggles_4(I).facLen_for_20 = dat(drg20(6:10),[2,6]);
    
    % figures
    Figure_1 = figure(I);
    subplot(3,2,1)
    plot(facLen_act_0,subAct_0,facLen_for_0,subFor_0);
    title(['Subject ',num2str(n(I)),': 0 drg']);
    xlabel('Fascicle Length');
    ylabel('Sub-Max');
    legend('Sub-Act','Sub-Force');
    grid on;
    
    subplot(3,2,2)
    plot(facLen_act_5,subAct_5,facLen_for_5,subFor_5);
    title(['Subject ',num2str(n(I)),': 5 drg']);
    xlabel('Fascicle Length');
    ylabel('Sub-Max');
    legend('Sub-Act','Sub-Force');
    grid on;
    
    subplot(3,2,3)
    plot(facLen_act_10,subAct_10,facLen_for_10,subFor_10);
    title(['Subject ',num2str(n(I)),': 10 drg']);
    xlabel('Fascicle Length');
    ylabel('Sub-Max');
    legend('Sub-Act','Sub-Force');
    grid on;
    
    subplot(3,2,4)
    plot(facLen_act_15,subAct_15,facLen_for_15,subFor_15);
    title(['Subject ',num2str(n(I)),': 15 drg']);
    xlabel('Fascicle Length');
    ylabel('Sub-Max');
    legend('Sub-Act','Sub-Force');
    grid on;
    
    subplot(3,2,5)
    plot(facLen_act_20,subAct_20,facLen_for_20,subFor_20);
    title(['Subject ',num2str(n(I)),': 20 drg']);
    xlabel('Fascicle Length');
    ylabel('Sub-Max');
    legend('Sub-Act','Sub-Force');
    grid on;
    
    
    
end

%% Force-length equation from Otten (1987)
addpath E:\Laptop\School\Interships\Zmat\Library
for I = 1:length(n)
    % grab subject data
    file = ['BAD_n',num2str(n(I)),'k',num2str(k),'.csv']; % File name
    data = csvread(file,1,0); % subject data
    
    % Separate and manipulate data
    A_magForce = V2F(data(1:25,5));  % magnitude of force at sub-activation (N)
    F_magForce = V2F(data(26:50,5)); % magnitude of force at sub-force (N)

%     [max_A_F, I_A_F] = max(A_magForce); % max force at act
%     [max_F_F, I_F_F] = max(F_magForce); % max force at force
% 
%     A_NF = (A_magForce/max_A_F)*100;     % normilised force at sub-act
%     force_NF = (F_magForce/max_F_F)*100; % normilised force at sub-force
% 
    A_Lf = data(1:25, 6); % facical length at sub-act
    F_Lf = data(26:50,6); % facical length at sub-force
% 
%     [max_A_Lf, I_A_Lf] = max(A_Lf); % max facical length at act
%     [max_F_Lf, I_F_Lf] = max(F_Lf); % max facical length  at force
% 
%     A_N_Lf = (A_Lf/max_A_Lf)*100; % normilised facical length at sub-act
%     F_N_Lf = (F_Lf/max_F_Lf)*100; % normilised facical length at sub-force

    A_NF = A_magForce; % do not normalise force at percent activation
    force_NF = F_magForce; % do not normalise force at percent force
    
    A_N_Lf = A_Lf; % do not normalise faciacl length at percent activation
    F_N_Lf = F_Lf; % do not normalise faciacl length at percent force
    
    % orgonise by sub-percentages
    count = 1;
    for i = 1:5:21

        A_20(count) = A_NF(i); % normal force at percent activation  
        A_40(count) = A_NF(i+1);
        A_60(count) = A_NF(i+2);
        A_80(count) = A_NF(i+3);
        A_100(count) = A_NF(i+4);

        Wiggles_4(I).A_20 = A_20; % save in struct
        Wiggles_4(I).A_40 = A_40;
        Wiggles_4(I).A_60 = A_60;
        Wiggles_4(I).A_80 = A_80;
        Wiggles_4(I).A_100 = A_100;
        
        Lf_A_20(count) = A_N_Lf(i); % normal facical length at percent activation 
        Lf_A_40(count) = A_N_Lf(i+1);
        Lf_A_60(count) = A_N_Lf(i+2);
        Lf_A_80(count) = A_N_Lf(i+3);
        Lf_A_100(count) = A_N_Lf(i+4);
        
        Wiggles_4(I).Lf_A_20 = Lf_A_20; % save in struct
        Wiggles_4(I).Lf_A_40 = Lf_A_40;
        Wiggles_4(I).Lf_A_60 = Lf_A_60;
        Wiggles_4(I).Lf_A_80 = Lf_A_80;
        Wiggles_4(I).Lf_A_100 = Lf_A_100;

        F_20(count) = force_NF(i); % normal force at percent force
        F_40(count) = force_NF(i+1);
        F_60(count) = force_NF(i+2);
        F_80(count) = force_NF(i+3);
        F_100(count) = force_NF(i+4);
        
        Wiggles_4(I).F_20 = F_20; % save in struct
        Wiggles_4(I).F_40 = F_40;
        Wiggles_4(I).F_60 = F_60;
        Wiggles_4(I).F_80 = F_80;
        Wiggles_4(I).F_100 = F_100;

        Lf_F_20(count) = F_N_Lf(i); % normal facical length at percent force 
        Lf_F_40(count) = F_N_Lf(i+1);
        Lf_F_60(count) = F_N_Lf(i+2);
        Lf_F_80(count) = F_N_Lf(i+3);
        Lf_F_100(count) = F_N_Lf(i+4);

        Wiggles_4(I).Lf_F_20 = Lf_F_20; % save in struct
        Wiggles_4(I).Lf_F_40 = Lf_F_40;
        Wiggles_4(I).Lf_F_60 = Lf_F_60;
        Wiggles_4(I).Lf_F_80 = Lf_F_80;
        Wiggles_4(I).Lf_F_100 = Lf_F_100;
        
        count = count+1;
    end 
    clear i count
    
    %  Activation Data Interplation
    long = max(Lf_A_20); %   interpolate act data 20
    short = min(Lf_A_20);
    Lf_A_20_new = long:-0.01:short;

    force_interp_20_A = interp1(Lf_A_20, A_20, Lf_A_20_new);
    Lf_A_20_new = Lf_A_20_new';

    
    FL_t_20 = 1:0.01/short:long/short; %   Filter act data
    %[freq, amp] = Freq_Amp_class(FL_t_20, force_interp_20_A);
    FL_dt_20 = FL_t_20(2) - FL_t_20(1);
    A_20_filt = filtmat_class(FL_dt_20', 15, force_interp_20_A');

    long = max(Lf_A_40); %   interpolate act data 40
    short = min(Lf_A_40);
    Lf_A_40_new = long:-0.01:short;

    force_interp_40_A = interp1(Lf_A_40, A_40, Lf_A_40_new);
    Lf_A_40_new = Lf_A_40_new';

    
    FL_t_40 = 1:0.01/short:long/short; %   Filter act data
    %[freq, amp] = Freq_Amp_class(FL_t_40, force_interp_40_A);
    FL_dt_40 = FL_t_40(2) - FL_t_40(1);
    A_40_filt = filtmat_class(FL_dt_40', 15, force_interp_40_A');

    
    long = max(Lf_A_60); %   interpolate act data 60
    short = min(Lf_A_60);
    Lf_A_60_new = long:-0.01:short;

    force_interp_60_A = interp1(Lf_A_60, A_60, Lf_A_60_new);
    Lf_A_60_new = Lf_A_60_new';

    
    FL_t_60 = 1:0.01/short:long/short; %   Filter act data
    %[freq, amp] = Freq_Amp_class(FL_t_60, force_interp_60_A);
    FL_dt_60 = FL_t_60(2) - FL_t_60(1);
    A_60_filt = filtmat_class(FL_dt_60', 15, force_interp_60_A');

    
    long = max(Lf_A_80); %   interpolate act data 80
    short = min(Lf_A_80);
    Lf_A_80_new = long:-0.01:short;

    force_interp_80_A = interp1(Lf_A_80, A_80, Lf_A_80_new);
    Lf_A_80_new = Lf_A_80_new';

    
    FL_t_80 = 1:0.01/short:long/short; %   Filter act data
    %[freq, amp] = Freq_Amp_class(FL_t_80, force_interp_80_A);
    FL_dt_80 = FL_t_80(2) - FL_t_80(1);
    A_80_filt = filtmat_class(FL_dt_80', 15, force_interp_80_A');

    
    long = max(Lf_A_100); %   interpolate act data 100
    short = min(Lf_A_100);
    Lf_A_100_new = long:-0.01:short;

    force_interp_100_A = interp1(Lf_A_100, A_100, Lf_A_100_new);
    Lf_A_100_new = Lf_A_100_new';

    
    FL_t_100 = 1:0.01/short:long/short; %   Filter act data
    %[freq, amp] = Freq_Amp_class(FL_t_100, force_interp_100_A);
    FL_dt_100 = FL_t_100(2) - FL_t_100(1);
    A_100_filt = filtmat_class(FL_dt_100', 15, force_interp_100_A');


    %  Force Data Interplation
    
    long = max(Lf_F_20); %   interpolate force data 20
    short = min(Lf_F_20);
    Lf_F_20_new = long:-0.01:short;

    force_interp_20_F = interp1(Lf_F_20, F_20, Lf_F_20_new);
    Lf_F_20_new = Lf_F_20_new';

    
    
    FL_t_20 = 1:0.01/short:long/short; %   Filter force data
    %[freq, amp] = Freq_Amp_class(FL_t_20, force_interp_20_F);
    FL_dt_20 = FL_t_20(2) - FL_t_20(1);
    F_20_filt = filtmat_class(FL_dt_20', 15, force_interp_20_F');

    
    long = max(Lf_F_40); %   interpolate force data 40
    short = min(Lf_F_40);
    Lf_F_40_new = long:-0.01:short;

    force_interp_40_F = interp1(Lf_F_40, F_40, Lf_F_40_new);
    Lf_F_40_new = Lf_F_40_new';

    
    FL_t_40 = 1:0.01/short:long/short; %   Filter force data
    %[freq, amp] = Freq_Amp_class(FL_t_40, force_interp_40_F);
    FL_dt_40 = FL_t_40(2) - FL_t_40(1);
    F_40_filt = filtmat_class(FL_dt_40', 15, force_interp_40_F');

    
    long = max(Lf_F_60); %   interpolate force data 60
    short = min(Lf_F_60);
    Lf_F_60_new = long:-0.01:short;

    force_interp_60_F = interp1(Lf_F_60, F_60, Lf_F_60_new);
    Lf_F_60_new = Lf_F_60_new';

    
    FL_t_60 = 1:0.01/short:long/short; %   Filter force data
    %[freq, amp] = Freq_Amp_class(FL_t_60, force_interp_60_F);
    FL_dt_60 = FL_t_60(2) - FL_t_60(1);
    F_60_filt = filtmat_class(FL_dt_60', 15, force_interp_60_F');

    
    long = max(Lf_F_80); %   interpolate force data 80
    short = min(Lf_F_80);
    Lf_F_80_new = long:-0.01:short;

    force_interp_80_F = interp1(Lf_F_80, A_80, Lf_F_80_new);
    Lf_F_80_new = Lf_F_80_new';

    
    FL_t_80 = 1:0.01/short:long/short; %   Filter force data
    %[freq, amp] = Freq_Amp_class(FL_t_80, force_interp_80_F);
    FL_dt_80 = FL_t_80(2) - FL_t_80(1);
    F_80_filt = filtmat_class(FL_dt_80', 15, force_interp_80_F');

    
    long = max(Lf_F_100); %   interpolate force data 100
    short = min(Lf_F_100);
    Lf_F_100_new = long:-0.01:short;

    force_interp_100_F = interp1(Lf_F_100, F_100, Lf_F_100_new);
    Lf_F_100_new = Lf_F_100_new';

    
    FL_t_100 = 1:0.01/short:long/short; %   Filter force data
    %[freq, amp] = Freq_Amp_class(FL_t_100, force_interp_100_F);
    FL_dt_100 = FL_t_100(2) - FL_t_100(1);
    F_100_filt = filtmat_class(FL_dt_100', 15, force_interp_100_F');

    %  Determine force-length properties
    [Fmmax_A_20, Lfopt_A_20, w_A_20, s_A_20, r_A_20] = FL_Fit_John(Lf_A_20_new', A_20_filt');
    [Fmmax_A_40, Lfopt_A_40, w_A_40, s_A_40, r_A_40] = FL_Fit_John(Lf_A_40_new', A_40_filt');
    [Fmmax_A_60, Lfopt_A_60, w_A_60, s_A_60, r_A_60] = FL_Fit_John(Lf_A_60_new', A_60_filt');
    [Fmmax_A_80, Lfopt_A_80, w_A_80, s_A_80, r_A_80] = FL_Fit_John(Lf_A_80_new', A_80_filt');
    [Fmmax_A_100, Lfopt_A_100, w_A_100, s_A_100, r_A_100] = FL_Fit_John(Lf_A_100_new', A_100_filt');

    [Fmmax_F_20, Lfopt_F_20, w_F_20, s_F_20, r_F_20] = FL_Fit_John(Lf_F_20_new', F_20_filt');
    [Fmmax_F_40, Lfopt_F_40, w_F_40, s_F_40, r_F_40] = FL_Fit_John(Lf_F_40_new', F_40_filt');
    [Fmmax_F_60, Lfopt_F_60, w_F_60, s_F_60, r_F_60] = FL_Fit_John(Lf_F_60_new', F_60_filt');
    [Fmmax_F_80, Lfopt_F_80, w_F_80, s_F_80, r_F_80] = FL_Fit_John(Lf_F_80_new', F_80_filt');
    [Fmmax_F_100, Lfopt_F_100, w_F_100, s_F_100, r_F_100] = FL_Fit_John(Lf_F_100_new', F_100_filt');

    G_Lfopt_A = Lfopt_A_100;
    G_Lfopt_F = Lfopt_F_100;

    % save Lfop data
    Wiggles_4(I).Lfopt_A_100 = Lfopt_A_100;
    Wiggles_4(I).Lfopt_A_80 = Lfopt_A_80;
    Wiggles_4(I).Lfopt_A_60 = Lfopt_A_60;
    Wiggles_4(I).Lfopt_A_40 = Lfopt_A_40;
    Wiggles_4(I).Lfopt_A_20 = Lfopt_A_20;
    
    Wiggles_4(I).Lfopt_F_100 = Lfopt_F_100;
    Wiggles_4(I).Lfopt_F_80 = Lfopt_F_80;
    Wiggles_4(I).Lfopt_F_60 = Lfopt_F_60;
    Wiggles_4(I).Lfopt_F_40 = Lfopt_F_40;
    Wiggles_4(I).Lfopt_F_20 = Lfopt_F_20;
    
    
    %   Activation data
    lf_A_20 = 0.5*Lfopt_A_20:0.1:1.5*Lfopt_A_20;
    FL_A_20 = forclen(lf_A_20, Lfopt_A_20, w_A_20, s_A_20, r_A_20);
    Lf_A_20_N = (Lf_A_20/G_Lfopt_A)*100;
    
%     figure(9+I)
%     plot(lf_A_20, FL_A_20*Fmmax_A_20, Lf_A_20_new, A_20_filt, 'g')

    lf_A_40 = 0.5*Lfopt_A_40:0.1:1.5*Lfopt_A_40;
    FL_A_40 = forclen(lf_A_40, Lfopt_A_40, w_A_40, s_A_40, r_A_40);
    Lf_A_40_N = (Lf_A_40/G_Lfopt_A)*100;
%     figure(9+I)
%     plot(lf_A_40, FL_A_40*Fmmax_A_40, Lf_A_40_new, A_40_filt, 'g')

    lf_A_60 = 0.5*Lfopt_A_60:0.1:1.5*Lfopt_A_60;
    FL_A_60 = forclen(lf_A_60, Lfopt_A_60, w_A_60, s_A_60, r_A_60);
    Lf_A_60_N = (Lf_A_60/G_Lfopt_A)*100;
%     figure(9+I)
%     plot(lf_A_60, FL_A_60*Fmmax_A_60, Lf_A_60_new, A_60_filt, 'g')

    lf_A_80 = 0.5*Lfopt_A_80:0.1:1.5*Lfopt_A_80;
    FL_A_80 = forclen(lf_A_80, Lfopt_A_80, w_A_80, s_A_80, r_A_80);
    Lf_A_80_N = (Lf_A_80/G_Lfopt_A)*100;
%     figure(9+I)
%     plot(lf_A_80, FL_A_80*Fmmax_A_80, Lf_A_80_new, A_80_filt, 'g')

    lf_A_100 = 0.5*Lfopt_A_100:0.1:1.5*Lfopt_A_100;
    FL_A_100 = forclen(lf_A_100, Lfopt_A_100, w_A_100, s_A_100, r_A_100);
    Lf_A_100_N = (Lf_A_100/G_Lfopt_A)*100;
%     figure(9+I)
%     plot(lf_A_100, FL_A_100*Fmmax_A_100, Lf_A_100_new, A_100_filt, 'g')

% Normilise Force and Facical lengths
    % Activation
    norm_lf_A_20 = (lf_A_20/Lfopt_A_100)*100;
    norm_FL_A_20 = ( (FL_A_20*Fmmax_A_20)/Fmmax_A_100 )*100;

    norm_lf_A_40 = (lf_A_40/Lfopt_A_100)*100;
    norm_FL_A_40 = ( (FL_A_40*Fmmax_A_40)/Fmmax_A_100 )*100;

    norm_lf_A_60 = (lf_A_60/Lfopt_A_100)*100;
    norm_FL_A_60 = ( (FL_A_60*Fmmax_A_60)/Fmmax_A_100 )*100;

    norm_lf_A_80 = (lf_A_80/Lfopt_A_100)*100;
    norm_FL_A_80 = ( (FL_A_80*Fmmax_A_80)/Fmmax_A_100 )*100;

    norm_lf_A_100 = (lf_A_100/Lfopt_A_100)*100;
    norm_FL_A_100 = ( (FL_A_100*Fmmax_A_100)/Fmmax_A_100 )*100;

    % save normal Lfop data
    Wiggles_4(I).norm_lf_A_20 = norm_lf_A_20;
    Wiggles_4(I).norm_lf_A_40 = norm_lf_A_40;
    Wiggles_4(I).norm_lf_A_60 = norm_lf_A_60;
    Wiggles_4(I).norm_lf_A_80 = norm_lf_A_80;
    Wiggles_4(I).norm_lf_A_100 = norm_lf_A_100;
    Wiggles_4(I).norm_FL_A_20 = norm_FL_A_20;
    Wiggles_4(I).norm_FL_A_40 = norm_FL_A_40;
    Wiggles_4(I).norm_FL_A_60 = norm_FL_A_60;
    Wiggles_4(I).norm_FL_A_80 = norm_FL_A_80;
    Wiggles_4(I).norm_FL_A_100 = norm_FL_A_100;
    
%     Figure_2 = figure(9+I);
%     plot(norm_lf_A_20, norm_FL_A_20, norm_lf_A_40, norm_FL_A_40,...
%     norm_lf_A_60, norm_FL_A_60, norm_lf_A_80, norm_FL_A_80, norm_lf_A_100,...
%     norm_FL_A_100);
%     title(['Subject ',num2str(n(I)),': Sub-max-act force length curves']);
%     ylabel('FDI Force, % of maximum at 100% activaiton');
%     xlabel('Fascicle Length, % of optimal');
%     legend('20%','40%', '60%', '80%', '100%');
%     grid on;

    %   Force data
    lf_F_20 = 0.5*Lfopt_F_20:0.1:1.5*Lfopt_F_20;
    FL_F_20 = forclen(lf_F_20, Lfopt_F_20, w_F_20, s_F_20, r_F_20);
    Lf_F_20_N = (Lf_F_20/G_Lfopt_F)*100;
%     figure(9+I)
%     plot(lf_F_20, FL_F_20*Fmmax_F_20, Lf_F_20_new, F_20_filt, 'g')

    lf_F_40 = 0.5*Lfopt_F_40:0.1:1.5*Lfopt_F_40;
    FL_F_40 = forclen(lf_F_40, Lfopt_F_40, w_F_40, s_F_40, r_F_40);
    Lf_F_40_N = (Lf_F_40/G_Lfopt_F)*100;
%     figure(9+I)
%     plot(lf_F_40, FL_F_40*Fmmax_F_40, Lf_F_40_new, F_40_filt, 'g')

    lf_F_60 = 0.5*Lfopt_F_60:0.1:1.5*Lfopt_F_60;
    FL_F_60 = forclen(lf_F_60, Lfopt_F_60, w_F_60, s_F_60, r_F_60);
    Lf_F_60_N = (Lf_F_60/G_Lfopt_F)*100;
%     figure(9+I)
%     plot(lf_F_60, FL_F_60*Fmmax_F_60, Lf_F_60_new, F_60_filt, 'g')

    lf_F_80 = 0.5*Lfopt_F_80:0.1:1.5*Lfopt_F_80;
    FL_F_80 = forclen(lf_F_80, Lfopt_F_80, w_F_80, s_F_80, r_F_80);
    Lf_F_80_N = (Lf_F_80/G_Lfopt_F)*100;
%     figure(9+I)
%     plot(lf_F_80, FL_F_80*Fmmax_F_80, Lf_F_80_new, F_80_filt, 'g')

    lf_F_100 = 0.5*Lfopt_F_100:0.1:1.5*Lfopt_F_100;
    FL_F_100 = forclen(lf_F_100, Lfopt_F_100, w_F_100, s_F_100, r_F_100);
    Lf_F_100_N = (Lf_F_100/G_Lfopt_F)*100;
%     figure(9+I)
%     plot(lf_F_100, FL_F_100*Fmmax_F_100, Lf_F_100_new, F_100_filt, 'g')
  
% Normilise Force and Facical lengths
    % Force
    norm_lf_F_20 = (lf_F_20/Lfopt_F_100)*100;
    norm_FL_F_20 = ( (FL_F_20*Fmmax_F_20)/Fmmax_F_100 )*100;

    norm_lf_F_40 = (lf_F_40/Lfopt_F_100)*100;
    norm_FL_F_40 = ( (FL_F_40*Fmmax_F_40)/Fmmax_F_100 )*100;

    norm_lf_F_60 = (lf_F_60/Lfopt_F_100)*100;
    norm_FL_F_60 = ( (FL_F_60*Fmmax_F_60)/Fmmax_F_100 )*100;

    norm_lf_F_80 = (lf_F_80/Lfopt_F_100)*100;
    norm_FL_F_80 = ( (FL_F_80*Fmmax_F_80)/Fmmax_F_100 )*100;

    norm_lf_F_100 = (lf_F_100/Lfopt_F_100)*100;
    norm_FL_F_100 = ( (FL_F_100*Fmmax_F_100)/Fmmax_F_100 )*100;
    
    % save normal Lfop data
    Wiggles_4(I).norm_lf_F_20 = norm_lf_F_20;
    Wiggles_4(I).norm_lf_F_40 = norm_lf_F_40;
    Wiggles_4(I).norm_lf_F_60 = norm_lf_F_60;
    Wiggles_4(I).norm_lf_F_80 = norm_lf_F_80;
    Wiggles_4(I).norm_lf_F_100 = norm_lf_F_100;
    Wiggles_4(I).norm_FL_F_20 = norm_FL_F_20;
    Wiggles_4(I).norm_FL_F_40 = norm_FL_F_40;
    Wiggles_4(I).norm_FL_F_60 = norm_FL_F_60;
    Wiggles_4(I).norm_FL_F_80 = norm_FL_F_80;
    Wiggles_4(I).norm_FL_F_100 = norm_FL_F_100;
    
%     Figure_3 = figure(100+I);
%     plot(norm_lf_F_20, norm_FL_F_20, norm_lf_F_40, norm_FL_F_40, norm_lf_F_60,...
%     norm_FL_F_60, norm_lf_F_80, norm_FL_F_80, norm_lf_F_100,...
%     norm_FL_F_100);
%     title(['Subject ',num2str(n(I)),': Sub-max-force force length curves']);
%     ylabel('FDI Force, % of maximum');
%     xlabel('Fascicle Length, % of optimal');
%     legend('20%','40%', '60%', '80%', '100%');
%     grid on;

   
end 
    
%% Save game

save('Wiggles_4.mat','Wiggles_4');
    
   









