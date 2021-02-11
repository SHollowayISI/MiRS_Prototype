%% MiRS - Simulation Setup
%{

    Sean Holloway
    MiRS Simulation Init File
    
    This file specifies simulation parameters for MiRS simulation.
    
%}

%% Simulation Parameter Setup

save_format.list = {'.png'};

% Radar simulation and processing setup
scenario.simsetup = struct( ...
    ...
    ... % Transceiver Trajectory Properties
    'radar_pos',    [0; 0; 0], ...              % Position of radar unit
    'radar_vel',    [0; 0; 0], ...              % Velocity of radar unit
    'radar_dir',    [0; 0], ...                 % Azimuth-Elevation normal of radar unit
    ...
    ... % Simulation Properties
    'num_frames',   1, ...                      % Number of radar frames to simulate
    'readout',      true, ...                   % Read out target data T/F
    ...
    ... % Processing Properties
    'par_test',     false, ...                  % Parallelize test loop
    'par_cfar',     true, ...                   % Parallelize CFAR detection
    'par_sim',      true, ...                   % Parallelize simulation
    'close_pool',   false, ...                  % Exit parallel pool at end
    ...
    'clear_cube',   false, ...
    'send_alert',   false, ...                  % Send email alert T/F
    'attach_zip',   false, ...
    'alert_address', 'sholloway@intellisenseinc.com', ...
    ...                                         % Email address for status updates
    'filename',     'Test_MiRS', ...            % Filename to save data as
    'save_format',  save_format, ...            % File types to save figures
    'save_figs',    false, ...                  % Save figures T/F
    'save_mat',     false, ...                  % Save mat file T/F
    'reduce_mat',   false);                     % Reduce mat file for saving


%% Start Parallel Pool

if (scenario.simsetup.par_cfar || scenario.simsetup.par_test)
    if(isempty(gcp('nocreate')))
        parpool;
    end
end

%% Test Mode Flag

test_mode = false;



