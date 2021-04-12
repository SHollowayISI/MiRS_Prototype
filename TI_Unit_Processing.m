
%% Housekeeping

close all;
clear variables;


%% Obtain input files

filepath = 'Input Data/Test 021621';
fold = dir(filepath);
k = 1;
files = {};
for i = 1:length(fold)
    if not(fold(i).isdir)
        if strcmp(fold(i).name((end-3):end), '.bin')
            files{k} = fold(i).name(1:end-4);
            k = k+1;
        end
    end
end

% names = {'Reflector Off', 'Active Reflector', 'Corner Reflector'};
% names = {'Test First Run', 'Test Second Run'};
% names = {'test'};
% names = {'Test', 'Far Reference', 'Far 1', 'Far 2', 'Far 3', 'Far 4', 'Far 5', ...
%     'Far 6', 'Far 7', 'Far 8', 'Far 9', 'Car', 'Far 10', 'Far Corner', 'Victor Approaching 1', ...
%     'Victor Approaching 2', 'Close Reference', 'Close 1', 'Close 2', 'Close 3', 'Close 4', ...
%     'Close 5', 'Close 6', 'Close 7', 'Close 8', 'Close 9', 'Close 10', 'Close Corner and Active', ...
%     'Close Corner'};
names = {'Far 1', 'Far 2', 'Far 3', 'Far 4', 'Far 5', ...
    'Far 6', 'Far 7', 'Far 8', 'Far 9', 'Far 10', 'Close 1', 'Close 2', 'Close 3', 'Close 4', ...
    'Close 5', 'Close 6', 'Close 7', 'Close 8', 'Close 9', 'Close 10'};

centerOfMass = {};
meanPosition = {};

for file_loop = 1:length(files)
    
    %% Settings
    
    f_c = 78e9;
    bw = 2e9;
    f_s = 20e6;
    t_ch = 102.4e-6;
    t_idle = 5e-6;
    
    c = physconst('Lightspeed');
    lambda = c / f_c;
    
    num_lanes = 2;
    num_samples = 4096;
    num_tx = 2;
    num_rx = 2;
    num_chirps = 232*2;
    num_frames = 10;
    complex_sampling = false;
    
    total_samps = num_samples*num_tx*num_rx*num_chirps*num_frames;
    if complex_sampling
        total_samps = total_samps * 2;
    end
    
    %% Parsing
    
    fid = fopen([filepath, '/', files{file_loop}, '.bin'], 'r');
    data = fread(fid, Inf, 'uint16=>single');
    fclose(fid);
    
    data = data - (data >= 2.^15) .* 2.^16;
    if complex_sampling
        data = reshape(data(1:total_samps), [8, total_samps/8]);
        data = permute(data(1:num_lanes,:) + 1i*data(4 + (1:num_lanes),:), [2 3 1]);
    else
        data = reshape(data(1:total_samps), [num_lanes, total_samps/num_lanes]);
        data = permute(data, [2 3 1]);
    end
    
    data = reshape(data, num_samples, [], num_rx);
    data_shaped = zeros(num_samples, num_chirps, num_tx*num_rx, num_frames);
    
    for r = 1:num_rx
        
        data_single_rx = reshape(data(:,:,r), num_samples, num_tx, num_chirps*num_frames);
        
        for t = 1:num_tx
            
            ind = (t-1)*num_rx + r;
            
            for fr = 1:num_frames
                
                fr_ind = (((fr - 1) * num_chirps) + 1):(fr * num_chirps);
                data_shaped(:,:,ind,fr) = squeeze(data_single_rx(:,t,fr_ind));
            end
        end
    end
    
    
    %% Radar cube processing
    
    N_r = 2^ceil(log2(num_samples));
    if complex_sampling
        range_cube = fft(hanning(size(data_shaped, 1)) .* data_shaped, N_r, 1);
        range_cube(end+1,:,:,:) = range_cube(1,:,:,:);
        range_res = (c/(2*bw));
        range_axis = ((0:N_r)-1)*range_res;
    else
        range_cube = fft(hanning(num_samples) .* data_shaped, N_r, 1);
        range_cube = range_cube(1:ceil(end/2),:,:,:);
        range_res = (c/(2*bw));
        range_axis = ((1:(N_r/2))-1)*range_res;
    end
    
    
    N_d = 2^ceil(log2(num_chirps));
    rd_cube = fftshift(fft(hanning(num_chirps)' .* range_cube, N_d, 2), 2);
%     rd_cube = fftshift(fft(range_cube, N_d, 2), 2);
    rd_cube(:,end+1,:,:) = rd_cube(:,1,:,:);
    
    vel_res = lambda/(2*(t_ch + t_idle)*num_chirps*num_tx);
    vel_axis = ((-N_d/2):(N_d/2))*vel_res;
    
    N_a = num_tx*num_rx;
    angle_cube = fftshift(fft(rd_cube, N_a, 3), 3);
    angle_cube(:,:,end+1,:) = angle_cube(:,:,1,:);
    
    pow_cube = abs(angle_cube).^2;
    static_cube = squeeze(pow_cube(:, ceil(end/2), ceil(end/2), :));
    
    %% Data Processing
    
    % CFAR Detection
    cfar = phased.CFARDetector( ...
        'Method',                   'CA', ...
        'NumGuardCells',            4, ...
        'NumTrainingCells',         8, ...
        'ProbabilityFalseAlarm',    1e-3, ...
        'ThresholdFactor',          'Auto', ...
        'OutputFormat',             'CUT result');
    
    if file_loop < 11
        ind = 1708:1736;
    else
        ind = 1042:1082;
    end
    
    int_axis = range_axis' .* static_cube;
    detect_list = cfar(static_cube, ind);
    detect_list = [zeros((ind(1)-1), num_frames); detect_list; zeros(size(static_cube,1) - ind(end), num_frames)];
    centerOfMass{file_loop} = nan(num_frames, 1);
    for fr = 1:num_frames
        centerOfMass{file_loop}(fr) = sum(int_axis(detect_list(:,fr) == 1)) / sum(static_cube(detect_list(:,fr) == 1));
    end
    meanPosition{file_loop} = mean(centerOfMass{file_loop}, 'omitnan');
    
    
    %% Plot results
    
   
%     figure('Name', 'Range Doppler Heat Map');
%     imagesc(vel_axis, range_axis, 10*log10(pow_cube(:,:,ceil(end/2), 1)));
%     set(gca, 'YDir', 'normal')
%     xlabel('Doppler Velocity [m/s]', 'FontWeight', 'bold')
%     ylim([0 20])
%     ylabel('Range [m]', 'FontWeight', 'bold');
    
    figure('Name', 'Zero Doppler Range Plot');
    for fr = 1:num_frames
        plot(range_axis, 10*log10(static_cube(:,fr)));
        hold on;
    end
    grid on;
    if file_loop < 11
        xlim([120 140])
    else
        xlim([70 90])
    end
    xlabel('Range [m]', 'FontWeight', 'bold');
    ylim([60 160])
    ylabel('Doppler Velocity [m/s]', 'FontWeight', 'bold')
    
%     figure('Name', 'Zero Doppler Range Plot Extended');
%     plot(range_axis, 10*log10(static_cube(:,1)));
%     grid on;
%     xlabel('Range [m]', 'FontWeight', 'bold');
%     xlim([0, range_axis(end)])
%     ylabel('Doppler Velocity [m/s]', 'FontWeight', 'bold')
    
    %% Save figures
    
    SaveFigures(names{file_loop}, 'Figures/Test 021621/Outside', '.png');
    close all;
    
    
    
end

%% All record processing

close all;

far_inds = 1:10;
close_inds = 11:20;

far_axis = (0:9)*0.1016;
close_axis = (0:9)*0.1016;

figure('Name', 'Far Position Range Estimation');
for n = 1:10
    scatter(far_axis(n)*ones(10,1), centerOfMass{far_inds(n)} - 128.244, '.', 'b');
    hold on;
end
scatter(far_axis, far_axis, 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Offset [m]');
xlim([-0.2, 1.2])
ylim([-0.2 1.2])

figure('Name', 'Far Position Average Range Estimation');
for n = 1:10
    scatter(far_axis(n), meanPosition{far_inds(n)} - 128.244, '.', 'b')
    hold on;
end
scatter(far_axis, far_axis, 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Offset [m]');
xlim([-0.2, 1.2])
ylim([-0.2 1.2])

figure('Name', 'Near Position Range Estimation');
for n = 1:10
    scatter(close_axis(n)*ones(10,1), centerOfMass{close_inds(n)} - 79.016, '.', 'b');
    hold on;
end
scatter(close_axis, close_axis, 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Offset [m]');
xlim([-0.2, 1.2])
ylim([-0.2 1.2])

figure('Name', 'Near Position Average Range Estimation');
for n = 1:10
    scatter(close_axis(n), meanPosition{close_inds(n)} - 79.016, '.', 'b')
    hold on;
end
scatter(close_axis, close_axis, 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Offset [m]');
xlim([-0.2, 1.2])
ylim([-0.2 1.2])

figure('Name', 'Far Position Range Error');
for n = 1:10
    scatter(far_axis(n)*ones(10,1), centerOfMass{far_inds(n)} - 128.244 - far_axis(n), '.', 'b');
    hold on;
%     scatter(far_axis(n), meanPosition{far_inds(n)} - 128.244 - far_axis(n), '+', 'b')
end
scatter(far_axis, zeros(size(far_axis)), 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Error [m]');
xlim([-0.2, 1.2])
ylim([-0.1 0.1])

figure('Name', 'Far Position Average Range Error');
for n = 1:10
    scatter(far_axis(n), meanPosition{far_inds(n)} - 128.244 - far_axis(n), '.', 'b')
    hold on;
end
scatter(far_axis, zeros(size(far_axis)), 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Error [m]');
xlim([-0.2, 1.2])
ylim([-0.1 0.1])

figure('Name', 'Near Position Range Error');
for n = 1:10
    scatter(close_axis(n)*ones(10,1), centerOfMass{close_inds(n)} - 79.016 - close_axis(n), '.', 'b');
    hold on;
end
scatter(close_axis, zeros(size(close_axis)), 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Error [m]');
xlim([-0.2, 1.2])
ylim([-0.2 0.2])

figure('Name', 'Near Position Average Range Error');
for n = 1:10
    scatter(close_axis(n), meanPosition{close_inds(n)} - 79.016 - close_axis(n), '.', 'b')
    hold on;
end
scatter(close_axis, zeros(size(close_axis)), 'o', 'r');
grid on;
xlabel('True Range Offset [m]');
ylabel('Measured Range Error [m]');
xlim([-0.2, 1.2])
ylim([-0.2 0.2])

SaveFigures('', 'Figures/Test 021621/Outside', '.png');
close all;

%% Function declarations

function [] = SaveFigures(save_name,fig_path, format)
%SAVEFIGURES Saves all open figures
%   Saves all open figures to fig_path directory, with name save_name plus
%   title of figure.

if ~exist(fig_path, 'dir')
    mkdir(fig_path)
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    saveas(FigHandle, fullfile(fig_path, [save_name, ' ', FigName, format]));
end

% Display update to command window
disp([format, ' Figures saved in ', save_name]);

end














