%% MiRS - Example Radar Initialization File
%{

    Sean Holloway
    MiRS Init File
    
    This file specifies radar parameters for MiRS simulation.
    
%}

%% Radar Parameter Setup

% Radar simulation and processing setup
scenario.radarsetup = struct( ...
    ...
    ... % Waveform Properties
    'f_c',          78.5e9, ...             % Operating frequency in Hz
    'f_s',          45e6, ...               % ADC sample frequency in Hz
    't_ch',         1.6e-3, ...             % Chirp duration in seconds
    'bw',           5.0e9, ...              % Chirp bandwidth in Hz
    'n_p',          64, ...                 % Number of (MIMO) chirps per CPI
    'drop_s',       7200, ...               % Number of samples to drop
    'cpi_fr',       5, ...                  % Number of CPI per frame
    ...
    ... % Antenna Array Properties
    'n_tx',         2, ...                  % Number of elements in horizontal Tx array
    'd_tx',         2, ...                  % Distance between Tx elements in wavelengths
    'n_rx',         4, ...                  % Number of elements in horizontal Rx array
    'd_rx',         0.5, ...                % Distance between Rx elements in wavelengths
    ...
    ... % Transceiver Properties
    'tx_pow',       db2pow(13 - 30), ...    % Transmit power per antenna in Watts
    'rf_sys_loss',  0, ...                  % RF system loss in dB
    'rx_nf',        12, ...                 % Rx noise figure in dB
    ...
    ... % Antenna Properties
    'ant_gain',     db2pow(10), ...         % Antenna gain in absolute
    'ant_cos_pow',  [1 10], ...             % Cosine power of antenna pattern
    ...
    ... % Processing Properties
    'r_win',        'hanning', ...          % Window for range processing
    'd_win',        'hanning', ...          % Window for doppler processing
    'az_win',       'none', ...          % Window for azimuth processing
    'n_az',         16, ...                 % Size of azimuth FFT
    'v_az_coeff',   -18.1246, ...           % Velocity-Bearing coupling in degree per meter per second
    ...
    ... % Detection Properties
    'detect_type',  'CFAR', ...             % Choose 'CFAR' or 'threshold'
    'thresh',       [], ...                 % Threshold in dB for threshold detection
    'CFAR_Pfa',     1e-6, ...             % CFAR false alarm probability
    'num_guard',    [3 1], ...              % Number of R-D guard cells for CFAR detection
    'num_train',    [15 2], ...              % Number of R-D training cells for CFAR detection
    'rm_group',     true, ...              % T/F remove closely grouped targets
    'rm_rad',       [1 Inf 5], ...          % Group radius in [m], [m/s], [deg]
    'dilate',       false, ...              % T/F dilate raw CFAR result to avoid duplicates 
    'dilate_shape', 'line', ...             % Shape to dilate raw CFAR results
    'dilate_args',  '3, 90', ...            % Arguments in strel function for dilation
    'det_m',        2);                     % M for m-of-n binary integration


%% Run Setup Scripts

% Set up Phased Array Toolbox system objects
scenario = PhasedSetup(scenario);





