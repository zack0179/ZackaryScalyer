%%Clean_Up%%%   Wiggling Fingers at the Force-Length Relationship: 
clear all   %   FDI, Architectural Influence on Force Potential.
clc         %   Zackary Scalyer
close all   %   The Pennsylvania State University - Berks
%%%%%%%%%%%%%   Spring_17: Dr. Ben Infantolino
%{
  Purpose - to analysis first dorsal interosseous (FDI) optimum fascicle 
  length at different sub-maximal force and sub-maximal percentages.  

  input
  -------
  ultrasound video of muscle fascicles (avi)
  EMG & force processed data (csv)

  output
  --------
  ActArch - Data set of FDI muscle architecture at sub-max activation
          - Column 1: Muscle thickness (Thickness)
          - Column 2: Pennation angle (Angle)
          - Column 3: Muscle fascicle length (Length)
  ForceArch - Data set of FDI muscle architecture at sub-max force
            - Column 1: Muscle thickness (Thickness)
            - Column 2: Pennation angle (Angle)
            - Column 3: Muscle fascicle length (Length)
  BAD - Single data set of all sub-max percentage for each joint angle
      - Column 1: Joint angle (Drg)
      - Column 2: Time (Time)
      - Column 3: Normalised EMG (NormEMG)
      - Column 4: Normalised Force (NormForce)
      - Column 5: Magnitude Force (MagForce)
      - Column 6: Muscle fascicle length (FacLength)
      - Column 7: Pennation angle (FacAngle)
      - Column 8: Muscle thickness (MusThick)

  Notes
  -------
  i. US mesuremets; Point 1 & 2, Thickness made in center of FDI belly,
                    Point 3 & 4, Supperficial apnorosis,
                    Point 5 & 6, Deep apnorosis
  ii. After reading data ensure only one value for each percentage. 
        - slect first occerence
        -7/17/2017 checks and corrections automated in filtering script
  ii. Fascicle length (Lfac) is not directly measured, it is calculated
      using pennation angle and muscle thickness according to Reeves and
      Narici (2003)


  Updates
  ---------
  i. 7-23-2017: Added varification of slected poitnts in cell 'FDI muscle
      architecture measurements'
  
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
%---------------------------%
% Enter investigator number %
%---------------------------%
N = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st MCP joint angles %%
%%%%%%%%%%%%%%%%%%%%%%%%%
MCP_drg = [0 5 10 15 20];
%%%%%%%%%%%%
% Percent %%
%%%%%%%%%%%%
Percent = 0.2;
disp('Finished Initialize variables');
%% Add path
%---------------------------%
% Enter subject folder path %
%---------------------------%
for i = 1:length(n)
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d',Drive, n(i)));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d\\n%d_out', Drive,n(i), n(i)));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d\\n%d_filtered', Drive,n(i), n(i)));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Research\\462w\\SubjectData\\n%d\\n%d_out\\investigator\\k%d', Drive,n(i), n(i), N));
    eval(sprintf('addpath %s:\\Laptop\\School\\Interships\\Zmat\\Library', Drive));
end
clear i

disp('Finished Add path');
%% Reading Data
%%%%%%%%%%%%%%%%%%%
% avi ultrasound %%
%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            fname_m = sprintf('S%d_ZS_%d.avi',n(I),MCP_drg(i));
            eval(sprintf('US_S%d_%d = VideoReader(fname_m);',n(I),MCP_drg(i)));    % US video
            eval(sprintf('num_US_S%d_%d = US_S%d_%d.NumberOfFrames;',...          % number of frames
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            eval(sprintf('FR_US_S%d_%d = US_S%d_%d.FrameRate;',...                % frame rate
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            eval(sprintf('dt_FR_US_S%d_%d = 1/FR_US_S%d_%d;',...                  % time interval
                n(I),MCP_drg(i),n(I),MCP_drg(i)));
            eval(sprintf('Time_US_S%d_%d = 0:dt_FR_US_S%d_%d:(num_US_S%d_%d*dt_FR_US_S%d_%d)-dt_FR_US_S%d_%d;',... % time array for each video
                n(I),MCP_drg(i),n(I),MCP_drg(i),n(I),MCP_drg(i),n(I),MCP_drg(i),n(I),MCP_drg(i)));       
    end
end
clear I i j fname_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove module level variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
        eval(sprintf('clear num_US_S%d_%d;', n,MCP_drg(i)));
        eval(sprintf('clear FR_US_S%d_%d;', n(I),MCP_drg(i)));
        eval(sprintf('clear dt_FR_US_S%d_%d;', n(I),MCP_drg(i)));
    end
end
clear I i j;
%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-Activation data %%
%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            fname = sprintf('SubAct_n%d_%ddrg.csv',n(I),MCP_drg(i));
            eval(sprintf('A_BAD_n%d_%ddrg = csvread(fname,1,0);',n(I),MCP_drg(i)));
    end
end
clear I i j fname;
%%%%%%%%%%%%%%%%%%%
% Sub-Force data %%
%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    for i = 1:length(MCP_drg)
            fname = sprintf('SubForce_n%d_%ddrg.csv',n(I),MCP_drg(i));
            eval(sprintf('F_BAD_n%d_%ddrg = csvread(fname,1,0);',n(I),MCP_drg(i)));
    end
end
clear I i j fname;

disp('Finished Reading Data');
%% FDI muscle architecture measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find US frame at sub-max activation occurrence %%
%          & measure FDI architecture            %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1:length(n)
    for i = 1:length(MCP_drg)
        Time_US_m = eval(sprintf('Time_US_S%d_%d',n(l),MCP_drg(i))); % grab US time
        SubAct_m = eval(sprintf('A_BAD_n%d_%ddrg',n(l),MCP_drg(i))); % grab EMG data
        K = zeros(size(SubAct_m,1),1);                               % initialise array for US index
        for ii = 1:length(SubAct_m) % Find index of US frame for each percentage of Sub-activation
            k = 0;                                  % index counter
            for jj = 1:length(Time_US_m) 
                if Time_US_m(jj) < SubAct_m(ii,1)   % if time of US frame is less than Sub-activation
                    k = k+1;                        % count
                else
                    K(ii) = k;                      % save US frame index
                end        
            end
        end
        I = 1;
        while I <= length(K) % for each US index 
            US_img_m = eval(sprintf('US_S%d_%d',n(l),MCP_drg(i))); % grab US frame
% play video up to the frame to be analysed
%             for J = 1:K(I)
%                 figure(1)
%                 pre_img_m = read(US_img_m,J); 
%                 imshow(pre_img_m,[])
%             end   
            img_m = read(US_img_m,K(I)); % read US frame
            figure(1)                    % figure to make architure mesuresments
            imshow(img_m)                % print US image
            [x,y] = ginput(6);           % slect six points(See Note i.)
            eval(sprintf('Athick_S%d_%d(I,1) = abs(y(2) - y(1))*(30/632);',...
                n(l),MCP_drg(i)));       % Point 1 & 2, Thickness made in center of FDI,
            eval(sprintf('ASupApon_S%d_%d(I,1) = atand(abs(y(4)-y(3))/abs(x(4)-x(3)));',...
                n(l),MCP_drg(i)));       % Point 3 & 4, Supperficial apnorosis,
            eval(sprintf('ADeepApon_S%d_%d(I,1) = atand(abs(y(6)-y(5))/abs(x(6)-x(5)));',...
                n(l),MCP_drg(i)));       % Point 5 & 6, Deep apnorosis
            eval(sprintf('APen_S%d_%d = ASupApon_S%d_%d + ADeepApon_S%d_%d;',... 
                n(l),MCP_drg(i),n(l),MCP_drg(i),n(l),MCP_drg(i))); % add subberficial and deep slopes
            eval(sprintf('AMusLen_S%d_%d = Athick_S%d_%d./sind(APen_S%d_%d);',... 
                n(l),MCP_drg(i),n(l),MCP_drg(i),n(l),MCP_drg(i))); % caculate angle of distal FDI, assumed pannation of facicals 
            lable_m = ['Investigator k',int2str(N),' MCP ',int2str(MCP_drg(i)),...
                'drg Act ',num2str(I*Percent)]; % lable for figure 2
            path_m = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                int2str(n(l)),'\n',int2str(n(l)),'_out\investigator\k',int2str(N),...
                '\','MCP_',int2str(MCP_drg(i)),'drg_',int2str(I),'Act','.fig']; % path for saving figure 2
            % Check Result
            figure(2)
            ax(1) = subplot(1,2,1); % axes one is US frame of FDI architure
            imshow(img_m)
            hold on
            h=plot(x(1:2),y(1:2),'-b',... % plot thickness and both aponeurosis
                x(3:4),y(3:4),'-r',... 
                x(5:6),y(5:6),'-r');
            set(h,'LineWidth',2);
            title(lable_m);
            ax(2) = subplot(1,2,2);
            img = imread('TryAgain.jpg'); % image to selcet to repet arch measures
            imshow(img); % plot imgage
            % a table of arch data
            eval(sprintf('Muscle_thick = Athick_S%d_%d(I,1);',n(l),MCP_drg(i)));
            eval(sprintf('Fac_length = AMusLen_S%d_%d(I,1);',n(l),MCP_drg(i)));
            t = table(Muscle_thick,Fac_length);
            uitable('Data',t{:,:},'ColumnName',t.Properties.VariableNames,...
                'RowName',t.Properties.RowNames,'Units', 'Normalized',...
                'Position',[.5 0 .272 .123]); % [left , boutum, width, height]
            [x] = ginput(1);
            axesSlect = gca;
            chooseData = (ax == axesSlect);
            if chooseData(1) == 1
                saveas(2,path_m); % save current mesures
                I = I+1; % itterate to next sub-activation
                close 1 2;
            else
                close 1 2;
            end  
        end
    end     
end

clear l i Time_US_m SubAct_m K ii k jj I US_img_m img_m pre_img_m lable_m path_m;
disp('finished sub-act');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find US frame at sub-max force occurrence %%
%        & measure FDI architecture         %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1:length(n)
    for i = 1:length(MCP_drg)
            Time_US_m = eval(sprintf('Time_US_S%d_%d',n(l),MCP_drg(i)));   % grab US time
            SubForce_m = eval(sprintf('F_BAD_n%d_%ddrg',n(l),MCP_drg(i))); % grab force data
            K = zeros(size(SubForce_m,1),1);                               % initialise array for US index
            for ii = 1:length(SubForce_m) % Find index of US frame for each percentage of Sub-force
                k = 0;                                  % index counter
                for jj = 1:length(Time_US_m) 
                    if Time_US_m(jj) < SubForce_m(ii,1) % if time of US frame is less than Sub-activation
                        k = k+1;                        % count
                    else
                        K(ii) = k;                      % save US frame index
                    end        
                end
            end
        I = 1;    
        while I <= length(K) % for each US index
            US_img_m = eval(sprintf('US_S%d_%d',n(l),MCP_drg(i))); % grab US frame
% play video up to the frame to be analysed
            for J = 1:K(I) 
               pre_img_m = read(US_img_m,J); 
               figure(1)
               imshow(pre_img_m,[])
               hold on;
            end   
            img_m = read(US_img_m,K(I)); % read US frame
            figure(1)                    % figure to make architure mesuresments
            imshow(img_m)                % print US image
            [x,y] = ginput(6);           % slect six points(See Note i.)
            eval(sprintf('Fthick_S%d_%d(I,1) = abs(y(2) - y(1))*(30/632);',...
                n(l),MCP_drg(i)));       % Point 1 & 2, Thickness made in center of FDI,
            eval(sprintf('FSupApon_S%d_%d(I,1) = atand(abs(y(4)-y(3))/abs(x(4)-x(3)));',...
                n(l),MCP_drg(i)));       % Point 3 & 4, Supperficial apnorosis,
            eval(sprintf('FDeepApon_S%d_%d(I,1) = atand(abs(y(6)-y(5))/abs(x(6)-x(5)));',...
                n(l),MCP_drg(i)));       % Point 5 & 6, Deep apnorosis
            eval(sprintf('FPen_S%d_%d = FSupApon_S%d_%d + FDeepApon_S%d_%d;',... 
                n(l),MCP_drg(i),n(l),MCP_drg(i),n(l),MCP_drg(i))); % add subberficial and deep slopes
            eval(sprintf('FMusLen_S%d_%d = Fthick_S%d_%d./sind(FPen_S%d_%d);',... 
                n(l),MCP_drg(i),n(l),MCP_drg(i),n(l),MCP_drg(i))); % caculate angle of distal FDI, assumed pannation of facicals 
            lable_m = ['Investigator k',int2str(N),' MCP ',int2str(MCP_drg(i)),...
                'drg Force ',num2str(I*Percent)]; % lable for figure 2
            path_m = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
                int2str(n(l)),'\n',int2str(n(l)),'_out\investigator\k',int2str(N),...
                '\','MCP_',int2str(MCP_drg(i)),'drg_',int2str(I),'Force','.fig']; % path for saving figure 2
            % Check Result
            figure(2)
            ax(1) = subplot(1,2,1); % axes one is US frame of FDI architure
            imshow(img_m)
            hold on
            h=plot(x(1:2),y(1:2),'-b',... % plot thickness and both aponeurosis
                x(3:4),y(3:4),'-r',... 
                x(5:6),y(5:6),'-r');
            set(h,'LineWidth',2);
            title(lable_m);
            ax(2) = subplot(1,2,2);
            img = imread('TryAgain.jpg'); % image to selcet to repet arch measures
            imshow(img); % plot imgage
            % a table of arch data
            eval(sprintf('Muscle_thick = Fthick_S%d_%d(I,1);',n(l),MCP_drg(i)));
            eval(sprintf('Fac_length = FMusLen_S%d_%d(I,1);',n(l),MCP_drg(i)));
            t = table(Muscle_thick,Fac_length);
            uitable('Data',t{:,:},'ColumnName',t.Properties.VariableNames,...
                'RowName',t.Properties.RowNames,'Units', 'Normalized',...
                'Position',[.5 0 .272 .123]); % [left , boutum, width, height]
            [x] = ginput(1);
            axesSlect = gca;
            chooseData = (ax == axesSlect);
            if chooseData(1) == 1
                saveas(2,path_m); % save current mesures
                I = I+1; % itterate to next sub-force
                close 1 2;
                clear Muscle_thick Fac_length
            else
                close 1 2;
                clear Muscle_thick Fac_length
            end  
        end
    end
end
clear l i Time_US_m SubForce_m K ii k jj I US_img_m img_m pre_img_m lable_m path_m
disp('finished sub-force');
%%%%%%%%%%%%%%%%%%%%%%
% Save Sub-max data %%
%%%%%%%%%%%%%%%%%%%%%%
for I = 1:length(n)
    k = 5;
    l=1;
    for i = 1:length(MCP_drg)
        eval(sprintf('Time(%d:%d,1) = A_BAD_n%d_%ddrg(:,1);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('NormEMG(%d:%d,1) = A_BAD_n%d_%ddrg(:,2);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('NormForce(%d:%d,1) = A_BAD_n%d_%ddrg(:,3);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('MagForce(%d:%d,1) = A_BAD_n%d_%ddrg(:,4);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('MusLength(%d:%d,1) = AMusLen_S%d_%d(:,1);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('FacAngle(%d:%d,1) = APen_S%d_%d(:,1);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('MusThick(%d:%d,1) = Athick_S%d_%d(:,1);',l,k,n(I),MCP_drg(i)));
        k=k+5;l=l+5;
    end
    for i = 1:length(MCP_drg)
        eval(sprintf('Time(%d:%d,1) = F_BAD_n%d_%ddrg(:,1);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('NormEMG(%d:%d,1) = F_BAD_n%d_%ddrg(:,2);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('NormForce(%d:%d,1) = F_BAD_n%d_%ddrg(:,3);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('MagForce(%d:%d,1) = F_BAD_n%d_%ddrg(:,4);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('MusLength(%d:%d,1) = FMusLen_S%d_%d(:,1);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('FacAngle(%d:%d,1) = Fthick_S%d_%d(:,1);',l,k,n(I),MCP_drg(i)));
        eval(sprintf('MusThick(%d:%d,1) = Fthick_S%d_%d(:,1);',l,k,n(I),MCP_drg(i)));
        k=k+5;l=l+5;
    end
    path = [Drive,':\Laptop\School\Interships\Zmat\Research\462w\SubjectData\n',...
        int2str(n(I)),'\n',int2str(n(I)),'_out\investigator\k',int2str(N),...
        '\BAD_n',int2str(n(I)),'k',int2str(N),'.csv'];
    Drg(1:25,1) = repelem(MCP_drg,5);
    Drg(26:50,1) = repelem(MCP_drg,5);
    t = table(Drg, Time, NormEMG, NormForce, MagForce, MusLength, FacAngle, MusThick);
    writetable(t,path);
end
clear I i l k path t Drg Time NormEMG NormForce MagForce MusLength FacAngle MusThick;

disp('Finished FDI muscle architecture measurements');
