files = uipickfiles;

silver_scaling = 65.534;
vel_scaling = 10;
knee_torque_scaling= 0.0276*178;
ankle_torque_scaling= 0.0276*200;

for i = 1:length(files)
    [fpath,fname,~] = fileparts(files{i});
    filevar = load(files{i});
    
    disp(fname);
    
    %Arrange in order of daq_data, locomode, no delay triggers, delay
    %triggers
    eval(['daq_data = filevar.',fname,'.daq.DAQ_DATA(:,[1:17 26:37]);']);
    eval(['loco_mode = filevar.',fname,'.pvd.SM_STATE_TAG(:,3) - 1005;']);
    eval(['combined_trig = filevar.',fname,'.pvd.SM_EVENT_LABEL;']);  
    
    daq_data(:,1:17) = daq_data(:,1:17) * silver_scaling;
    daq_data(:,[2 5]) = daq_data(:,[2 5]) * vel_scaling;
    daq_data(:,[3 15]) = daq_data(:,[3 15]) * knee_torque_scaling;
    daq_data(:,[6 16]) = daq_data(:,[6 16]) * ankle_torque_scaling;
       
    nondelay_trig = combined_trig;
    delay_trig = combined_trig;
    
    for j = 1:length(combined_trig)
        if j < length(combined_trig) - 3
            if combined_trig(j) > 10000
                nondelay_trig(j-3) = combined_trig(j)/10;
                nondelay_trig(j) = 0;
            end
        end
    end
    nondelay_trig(end-2:end) = 0;
    
    for j = 1:length(combined_trig)
        if j > 3
            if combined_trig(j) > 10000
                delay_trig(j-3) = 0;
            end
        end
    end
    todelay_trig = find(delay_trig < 10000 & delay_trig > 0);
    todelay_trig = todelay_trig(todelay_trig < length(combined_trig)-3);
    delay_trig(todelay_trig+3) = delay_trig(todelay_trig)*10;
    delay_trig(find(delay_trig < 10000)) = 0;
    
    loco_mode_daq = [];
    nondelay_trig_daq = [];
    delay_trig_daq = [];
    for j = 1:length(nondelay_trig)
        loco_mode_daq = [loco_mode_daq; repmat(loco_mode(j),30,1)];
        nondelay_trig_daq = [nondelay_trig_daq; repmat(nondelay_trig(j),30,1)];
        delay_trig_daq = [delay_trig_daq; repmat(delay_trig(j),30,1)];
    end
    
    sensor_startind = [26 32];
    
    % Put IMU_data in cells based on number of sensors and apply scaling
    for i = 1:length(sensor_startind)
        % Acc scaling: (UINT - 32768)/8192
        % Gyro scaling: (UINT - 32768)/65.536
        Accel = eval(['(double(filevar.',fname,'.daq.daqUINT16(:,sensor_startind(i):sensor_startind(i)+2)) - 32768)/8192;']);
        Gyro = eval(['(double(filevar.',fname,'.daq.daqUINT16(:,sensor_startind(i)+3:sensor_startind(i)+5)) - 32768)/65.536;']);
        IMU_data{i} = [Accel Gyro];
    end
    
    Fs = 1000;
    alpha = 0.99;
    NetAcc_low = 0.98;
    NetAcc_high = 1.02;
    
    for i = 1:length(IMU_data)
        Ax = IMU_data{i}(:,1);
        Ay = IMU_data{i}(:,2);
        Az = IMU_data{i}(:,3);
        Gy = IMU_data{i}(:,4);
        Gz = IMU_data{i}(:,5);
        Gx = IMU_data{i}(:,6);
        
        NetAcc = sqrt(Ax.^2 + Az.^2);
        
        pitch = atan2d(Az,Ax.^2+Ay.^2);
        
        % Define initial pose
        unfilt_angle(1) = pitch(1);
        filt_angle(1) = pitch(1);
        
        for j = 2:size(Gy,1)
            % Apply complementary filter when accel. reading is in normal range
            unfilt_angle(j) = unfilt_angle(j-1) + Gy(j)*1/Fs;
            if NetAcc(j) > NetAcc_low && NetAcc(j) < NetAcc_high
                filt_angle(j) = alpha*(filt_angle(j-1) + Gy(j)*1/Fs) + (1-alpha)*pitch(j);
            else
                filt_angle(j) = filt_angle(j-1) + Gy(j)*1/Fs;
            end
        end
        
        % Subtract initial pose
        unfilt_angle = unfilt_angle - pitch(1);
        filt_angle = filt_angle - pitch(1);
        pitch = pitch - pitch(1);
        
        IMU_data{i} = [IMU_data{i} pitch unfilt_angle' filt_angle'];
        
        pitch = [];
        unfilt_angle = [];
        filt_angle = [];
    end
    
    shank_angle = IMU_data{1}(:,9);
    thigh_angle = IMU_data{2}(:,9);
    
    python_input = [daq_data(:,1:17) IMU_data{1}(:,1:6) IMU_data{2}(:,1:6) shank_angle thigh_angle loco_mode_daq nondelay_trig_daq delay_trig_daq];
    chanlabels = {'Knee Angle','Knee Vel','Knee Current','Ankle Angle','Ankle Vel','Ankle Current','VU Ax','VU Ay','VU Az','VU Gx','VU Gz','VU Gy','Shank Angle','Thigh Angle','Knee Ref','Ankle Ref','Load','Shank Ax','Shank Ay','Shank Az','Shank Gy','Shank Gz','Shank Gx','Thigh Ax','Thigh Ay','Thigh Az','Thigh Gy','Thigh Gz','Thigh Gx','Contra Shank','Contra Thigh'};
    
    save(['Z:\Lab Member Folders\Blair Hu\Contralateral Prosthesis Control 2017\TF01_Goldie04_121917_Dev\DATA\FORPY\',fname,'_pydaq.mat'],'python_input','chanlabels');
end