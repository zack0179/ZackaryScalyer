%%Clean_Up%%%   Zacary Scalyer
clearvars   %   Dynamic test data
clc         %   PSU Berks Summer 2017
close all   %   Biomechanics Internship
%%%%%%%%%%%%%   Adviser - Dr. Joe Mahoney
%{
Purpose: To compare a cost effective Phidget force plate (PFP) with Bertec
    6090 and/or Bertec 4060 force plate (BFP6090 / BFP4060) during quiet 
    stance.
Proceger: The Phidget Force Plate (PFP) was placed on the Bertec 6090 Force 
    Plate (BFP6090) and recored three quiet stance trials with feet together  
    and hands on the illiac crest. PFP data was colected with a computer 
    using python 2.7 while the BFP6090 data was colected with the Vicon 
    Nexis II system. The two signals were syniced with a Ardunio
    differenctal trigger.
UpData: 
    07-19-2017: - Bertec6090 Adj COP, h was meaured with a digital digital 
                    caliper in mm, converted to m.
                - Phidget COP was in m, converted to mm.
                - Inverse Y axes seems to have allined the overlay of PFP
                    and Adj Bertec COP
    07-21-2017: - One loop for a single struct containing all subject data.
%}

%% Data collection parameters
clc;
%%%%%%%%%%%%%%%%%%%
% Subject Number %%
%%%%%%%%%%%%%%%%%%%
n = 20;

%%%%%%%%%%%%
% Subject %%
%%%%%%%%%%%%
Subject = [];
for i = 1:n
    if i < 10
        Subject = [Subject; strjoin({'AB0',num2str(i)},'')];
    else
        Subject = [Subject; strjoin({'AB',num2str(i)},'')];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% Names of condtions %%
%%%%%%%%%%%%%%%%%%%%%%%
Condition = {'Double_Open','Double_Closed','Single_Open', 'Tandem_Open'};

%%%%%%%%%%%%%%%
% Conditions %%
%%%%%%%%%%%%%%%
Test = [1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data interval of interest %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Start = 5;       % start point (seconds)
Stop = 30+Start; % end point

disp('Finished: Data collection parameters');
disp('Open: Dynamic folder in MATLAB Current Folder');
%% Calibration parameters
clc;

%%%%%%%%%%%%%%%%%%%%%%%%
% Loadcell caibration %%
%%%%%%%%%%%%%%%%%%%%%%%%
load('Calibration3.mat'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Possions of sesors on PFP %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensor 0 is the origin
r01 = [479.18e-3; 180-0]; % r [m] & theta [deg] from origin to test point 1
r02 = [479.88e-3,180-61]; % r [m] & theta [deg] from origin to test point 2
r01 = r01(1)*[cosd(r01(2)),sind(r01(2))]; % convert to cartesian in [m]
r02 = r02(1)*[cosd(r02(2)),sind(r02(2))]; % convert to cartesian in [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for x,y center of mass on PFP %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xpos = @(f0,f1,f2) (f1*r01(1)+f2*r02(1) )./(f0+f1+f2);
ypos = @(f0,f1,f2) (f1*r01(2)+f2*r02(2) )./(f0+f1+f2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for x,y center of mass on BFP %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for x,y of a force applied to the Bertec plates at a distence h

Xp = @(h,Fx,My,Fz) (-h.*Fx-My)./Fz;
Yp = @(h,Fy,Mx,Fz) (-h.*Fy+Mx)./Fz;

h = mean([109.43,112.55,108.78])/1000; % m form the top of the Bertec to the application of force on PFP

disp('Finished: Calibration perameters');
disp('Open subject folder in MATLAB Current Folder for next cell');

%% The Struct

cop_dat = struct('Subject',{}, 'Condition', {}, 'Test', [], 'Sync', [], ... 
    'PFP_Time', [],'PFP_LcForce_0', [], 'PFP_LcForce_1', [],'PFP_LcForce_2', [],...
    'PFP_NetForce', [], 'PFP_X_cordant', [],'PFP_Y_cordant', [],...
    'PFP_radii', [], 'PFP_PathLength', [],'PFP_eigenVal', [], ...                
    'PFP_eigenVec', [], 'PFP_AxesLengths', [], 'PFP_Area', [],...
    'BFP_Time', [],'BFP_Force_Z', [], 'BFP_X_raw', [], 'BFP_Y_raw', [],...
    'BFP_X_cordant', [],'BFP_Y_cordant', [], 'BFP_radii', [], ...
    'BFP_PathLength', [], 'BFP_eigenVal', [],'BFP_eigenVec', [], ...
    'BFP_AxesLengths', [], 'BFP_Area', [],'PFP_maxvel',[],'BFP_maxvel',[]);

disp('Finished Initslisesing Matrix of Structs');
%% Read analysis and compile data from files
close all; clc;
% List of subfolders in subject data
listDir = subdir('E:\Laptop\School\Interships\Zmat\Research\ForcePlate_Phidget\Aquire4\Dynamic\Subject_data');
count = 0;
missingData = 0;
for J = 1:length(listDir) 
    for I = 1:length(Condition)
        for i = 1:length(Test)
            PFP_path = strjoin([listDir(J),'\','Phidget_',Condition(I),'_',num2str(Test(i)),'.csv'],''); % path and file name for Phidet force plate
            BFP_path = strjoin([listDir(J),'\','Bertec6090_',Condition(I),'_',num2str(Test(i)),'.csv'],''); % path and file name for Bertec force plate
            try
                PFP_dat = csvread(PFP_path,1,0); % atempt to read in data
                BFP_dat = csvread(BFP_path,1,0);   

            catch
                missingData = missingData +1;
                continue % if file does not exist, itterate to next listDir(J)
            end
            % add subject and condtion data to struct
            count = count +1;
            cop_dat(count).Subject = Subject(J,:);
            cop_dat(count).Test = Test(i);
            cop_dat(count).Condition = Condition(I);

            % Phidget data
            inds0 = find(PFP_dat(:,1)==0); % sensor 0
            inds1 = find(PFP_dat(:,1)==1); % sensor 1
            inds2 = find(PFP_dat(:,1)==2); % sensor 2

            if ~isempty(inds0)
                dat0 = PFP_dat(inds0,2);
                tim0 = PFP_dat(inds0,3);
                tspan=0:1/120:min(tim0(end)); % uniform time for all sensors 
            end
            if ~isempty(inds1)
                dat1 = PFP_dat(inds1,2);
                tim1 = PFP_dat(inds1,3);
                tspan=0:1/120:min([tim0(end) tim1(end)]); % uniform time for all sensors 
            end
            if ~isempty(inds2)
                dat2 = PFP_dat(inds2,2);
                tim2 = PFP_dat(inds2,3);
                tspan=0:1/120:min([tim0(end) tim1(end) tim2(end)]); % uniform time for all sensors 
            end


            dat0trim = interp1(tim0,dat0,tspan);    % resample to common time 
            dat1trim = interp1(tim1,dat1,tspan);
            dat2trim = interp1(tim2,dat2,tspan);

            f0 = [dat0trim' ones(length(dat0trim),1)]*V0; % find force on each sensor
            f1 = [dat1trim' ones(length(dat1trim),1)]*V1;
            f2 = [dat2trim' ones(length(dat2trim),1)]*V2;

            NetForce = f0+f1+f2; % Force (N)

            mark = [];
            neton = 50;
            while isempty(mark)
                [mark] = find(NetForce >= neton + median(NetForce,'omitnan'),1); % mark point of synchronisation
                neton = neton - 1;
                if neton < 5
                    disp(BFP_path);disp(PFP_path);
                    disp('Sync Point Not found in Phidget data');
                    return
                end
            end

            % Check data
%             figure(1)
%             subplot(1,2,1);
%             plot(NetForce);
%             subplot(1,2,2);
%             plot(BFP_dat(:,4));
            
            [b,a] = butter(4,5/120); % 4th order transfer function coefficients, cut off/sample rate  

            try
                NetForce =  filtfilt(b,a,NetForce(mark+Start*120:mark+Stop*120)); % filter and trim interval of synchronised data
                f0 =  filtfilt(b,a,f0(mark+Start*120:mark+Stop*120));
                f1 =  filtfilt(b,a,f1(mark+Start*120:mark+Stop*120));
                f2 =  filtfilt(b,a,f2(mark+Start*120:mark+Stop*120));
                cop_dat(count).Sync = 1;
            catch
                NetForce =  filtfilt(b,a,NetForce(Start*120:Stop*120)); % filter and trim interval of synchronised data
                f0 =  filtfilt(b,a,f0(Start*120:Stop*120));
                f1 =  filtfilt(b,a,f1(Start*120:Stop*120));
                f2 =  filtfilt(b,a,f2(Start*120:Stop*120));
                disp(PFP_path); disp('Phidget 30s widow is outside data range');
                disp('Data is not allined');
                cop_dat(count).Sync = 0;
            end
                

            XP = xpos(f0,f1,f2)/1000; % convert three force measurements to xy cop [mm]
            YP = ypos(f0,f1,f2)/1000; 

            XPnorm = (XP - mean(XP))*1000; % Normalised xy cop [mm]
            YPnorm = (YP - mean(YP))*1000;

            radii = sqrt( (XPnorm.^2) + (YPnorm.^2 )); % radius of COP
            ds = sqrt(sum(diff([XPnorm,YPnorm]).^2,2)); % change in position 
            pathLen = sum(ds);
            cop_dat(count).PFP_maxvel = max(ds)*120; % max velocity in m/s

            % principal component analysis
            CXY = [XPnorm, YPnorm]; 
            Cp = (CXY'*CXY)./(length(CXY)-1); 
            [V, pc] = eig(Cp);   % V - eigenVectors, pc - eigenValues 
            pc = sort(diag(pc)); % eigenValues 

            % Ellipses fit
            conf = 0.95;                                    % confidence of ellipse
            [n,p] = size(CXY);                              % n-number of points, p- dimensions in the ellipse
            f95 = finv(0.95,p,n-p)*(n-1)*p*(n+1)/n/(n-p);   % 'F 95 percent point function'; F-CDF with inverse F-function?
            saxes = sort(sqrt(pc*f95),'descend');           % 'semi-axes lengths'
            area = pi^(p/2)/gamma(p/2+1)*prod(saxes);% Area of 0.95 ellipse fit

            % Phidet data to struct
            cop_dat(count).PFP_Time = tspan;
            cop_dat(count).PFP_LcForce_0 = f0;
            cop_dat(count).PFP_LcForce_1 = f1;
            cop_dat(count).PFP_LcForce_2 = f2;
            cop_dat(count).PFP_NetForce = NetForce;
            cop_dat(count).PFP_X_cordant = -XPnorm;
            cop_dat(count).PFP_Y_cordant = YPnorm;
            cop_dat(count).PFP_radii = radii;
            cop_dat(count).PFP_PathLength = pathLen;
            cop_dat(count).PFP_eigenVal = pc;
            cop_dat(count).PFP_eigenVec = V;
            cop_dat(count).PFP_AxesLengths = saxes;
            cop_dat(count).PFP_Area = area;

            % Bertec data

            mark = [];
            neton = 50;
            while isempty(mark)
                [mark] = find(abs(BFP_dat(:,4)) > neton + median(abs(BFP_dat(:,4))), 1); % find point of synchronisation
                neton = neton - 1;
                if neton < 10
                    disp(BFP_path);disp(PFP_path);
                    disp('Sync Point Not found in Bertec data');
                    return
                end
            end

            try
                trim_BFP_dat = BFP_dat(mark+Start*1000:mark+Stop*1000,:);  % truncate data
                if cop_dat(count).Sync ~= 0
                    cop_dat(count).Sync = 1;
                end
            catch
                trim_BFP_dat = BFP_dat(Start*1000:Stop*1000,:);
                disp(BFP_path); disp('Bertec 30s widow is outside data range');
                disp('Data is not allined');
                cop_dat(count).Sync = 0;
            end


            X_raw = trim_BFP_dat(:,8);        % X center of pressure [mm]
            Y_raw = trim_BFP_dat(:,9);        % Y center of pressure [mm]

            Fz = trim_BFP_dat(:,4); % Froce in Z
            Fx = trim_BFP_dat(:,2); % Froce in x
            Fy = trim_BFP_dat(:,3); % Froce in y
            Mx = trim_BFP_dat(:,5); % Moment in x
            My = trim_BFP_dat(:,6); % Moment in y

            XP = Xp(h,Fx,My,Fz); % x,y cordents for force at dist h 
            YP = Yp(h,Fy,Mx,Fz);

            % create time array
            time = 0:1/1000:(size(Fz,1)-1)/1000;    % Bertec time (seconds) 1000Hz
            t = 0:1/120:(size(Fz,1)-1)/1000;        % interpal time 120Hz

            % resample to common time 
            interp_XP = interp1(time, XP, t)';       % x cordent for center of pressure (COP)
            interp_YP = interp1(time, YP,  t)';      % Y cordent 
            interp_force = interp1(time, Fz, t)';    % force in Z-axis (N)
            interp_X_raw = interp1(time, X_raw, t)'; 
            interp_Y_raw = interp1(time, Y_raw, t)'; 

            [b,a] = butter(4,5/120); % 4th order transfer function coefficients, cut off/sample rate  

            filter_XP =  filtfilt(b,a,interp_XP); % filter and trim interval of synchronised data
            filter_YP =  filtfilt(b,a,interp_YP);
            filter_force =  filtfilt(b,a,interp_force);
            filter_X_raw =  filtfilt(b,a,interp_X_raw);
            filter_Y_raw =  filtfilt(b,a,interp_Y_raw);

            % mean zero raw cop
            X_raw = (filter_X_raw - mean(filter_X_raw));
            Y_raw = (filter_Y_raw - mean(filter_Y_raw));

            % mean zero adjed cop
            xpnorm = filter_XP - mean(filter_XP); 
            ypnorm = filter_YP - mean(filter_YP);  

            radii = sqrt( (xpnorm.^2) + (ypnorm.^2 )); % radius of COP
            ds = sqrt(sum(diff([xpnorm,ypnorm]).^2,2)); % change in position 
            pathLen = sum(ds);
            cop_dat(count).BFP_maxvel = max(ds)*120; % max velocity in m/s

            % principal component analysis
            CXY = [xpnorm, ypnorm]; 
            Cp = (CXY'*CXY)./(length(CXY)-1); 
            [V, pc] = eig(Cp); % V - eigenVectors, pc - eigenValues 
            pc = sort(diag(pc)); % eigenValues 

            % Ellipses fit
            conf = 0.95; % confidence of ellipse
            [n,p] = size(CXY); % n-number of points, p- dimensions in the ellipse
            f95 = finv(0.95,p,n-p)*(n-1)*p*(n+1)/n/(n-p); % 'F 95 percent point function'; F-CDF with inverse F-function?
            saxes = sort(sqrt(pc*f95),'descend'); % 'semi-axes lengths'
            area = pi^(p/2)/gamma(p/2+1)*prod(saxes); % volume of 0.95 hyper-ellipsoid fit

            % Bertec data to struct
            cop_dat(count).BFP_Time = t;
            cop_dat(count).BFP_Force_Z = filter_force;
            cop_dat(count).BFP_X_raw = X_raw;
            cop_dat(count).BFP_Y_raw = Y_raw;
            cop_dat(count).BFP_X_cordant = xpnorm;
            cop_dat(count).BFP_Y_cordant = ypnorm;
            cop_dat(count).BFP_radii = radii;
            cop_dat(count).BFP_PathLength = pathLen;
            cop_dat(count).BFP_eigenVal = pc;
            cop_dat(count).BFP_eigenVec = V;
            cop_dat(count).BFP_AxesLengths = saxes;
            cop_dat(count).BFP_Area = area;
        end
    end
end
clear i I j PFP_path PFP_dat inds0 inds1 inds2 dat0trim dat1trim dat2trim...
    dat0 dat1 dat2 f0 f1 f2 NetForce mark a b XP YP XPnorm YPnorm radii...
    tspan tim0 tim1 tim2 pathLen CXY Cp V pc conf n p f95 saxes area...
    i I j BFP_path BFP_dat mark X_raw Y_raw Fz Fx Fy Mx My XP YP mark time...
    t interp_XP interp_YP interp_force interp_X_raw interp_Y_raw xpnorm...
    ypnorm radii pathLen CXY Cp pc V conf n p f95 saxes area

disp('Finished: Read analysis and compile data from files');

%% Save game

%%%%%%%%%%%%%%%
% Structures %%
%%%%%%%%%%%%%%%
save('cop_data.mat','cop_dat');





%%

