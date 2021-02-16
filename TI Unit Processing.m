
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

names = {'Reflector Off', 'Active Reflector', 'Corner Reflector'};

for file_loop = 1:length(files)
    
    %% Settings
    
    f_c = 78e9;
    bw = 2e9;
    f_s = 40e6;
    t_ch = 102.4e-6;
    t_idle = 5e-6;
    
    c = physconst('Lightspeed');
    lambda = c / f_c;
    
    num_lanes = 2;
    num_samples = 4096;
    num_tx = 2;
    num_rx = 2;
    num_chirps = 255;
    num_frames = 1;
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
        data = reshape(data(1:total_samps), [2*num_lanes, total_samps/(2*num_lanes)]);
        data = permute(data(1:num_lanes,:) + 1i*data((num_lanes+1):(2*num_lanes),:), [2 3 1]);
    else
        data = reshape(data(1:total_samps), [num_lanes, total_samps/num_lanes]);
        data = permute(data, [2 3 1]);
    end
    
    data = reshape(data, num_samples, [], num_rx);
    data_shaped = zeros(num_samples, num_chirps*num_frames, num_tx*num_rx);
    
    for r = 1:num_rx
        
        data_single_rx = reshape(data(:,:,r), num_samples, num_tx, num_chirps*num_frames);
        
        for t = 1:num_tx
            ind = (t-1)*num_rx + r;
            
            data_shaped(:,:,ind) = squeeze(data_single_rx(:,t,:));
        end
    end
    
    
    %% Radar cube processing
    
    N_r = 2^ceil(log2(num_samples));
    range_cube = fft(hanning(num_samples) .* data_shaped, N_r, 1);
    if complex_sampling
        range_cube(end+1,:,:) = range_cube(1,:,:);
        range_res = (c/(2*bw));
        range_axis = ((0:N_r)-1)*range_res;
    else
        range_cube = range_cube(1:ceil(end/2),:,:);
        range_res = (c/(2*bw));
        range_axis = ((1:(N_r/2))-1)*2*range_res;
    end
    
    
    N_d = 2^ceil(log2(num_chirps*num_frames));
    rd_cube = fftshift(fft(hanning(num_chirps*num_frames)' .* range_cube, N_d, 2), 2);
    rd_cube(:,end+1,:) = rd_cube(:,1,:);
    
    vel_res = lambda/(2*(t_ch + t_idle)*num_chirps*num_tx);
    vel_axis = ((-N_d/2):(N_d/2))*vel_res;
    
    N_a = num_tx*num_rx;
    angle_cube = fftshift(fft(rd_cube, N_a, 3), 3);
    angle_cube(:,:,end+1) = angle_cube(:,:,1);
    
    pow_cube = abs(angle_cube).^2;
    
    
    %% Plot results
    
    figure('Name', 'Range Doppler Heat Map');
    imagesc(vel_axis, range_axis, 10*log10(pow_cube(:,:,ceil(end/2))));
    set(gca, 'YDir', 'normal')
    xlabel('Doppler Velocity [m/s]', 'FontWeight', 'bold')
    ylim([0 20])
    ylabel('Range [m]', 'FontWeight', 'bold');
    
    figure('Name', 'Zero Doppler Range Plot');
    plot(range_axis, 10*log10(pow_cube(:,ceil(end/2),ceil(end/2))));
    grid on;
    xlim([0 20])
    xlabel('Range [m]', 'FontWeight', 'bold');
    ylim([60 160])
    ylabel('Doppler Velocity [m/s]', 'FontWeight', 'bold')
    
    figure('Name', 'Zero Doppler Range Plot Extended');
    plot(range_axis, 10*log10(pow_cube(:,ceil(end/2),ceil(end/2))));
    grid on;
    xlabel('Range [m]', 'FontWeight', 'bold');
    xlim([0, range_axis(end)])
    ylabel('Doppler Velocity [m/s]', 'FontWeight', 'bold')
    
    %% Save figures
    
    SaveFigures(names{file_loop}, 'Figures/Test 021621', '.png');
    close all;
    
    
    
end


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














