function [cube] = SignalProcessing(scenario)
%SIGNALPROCESSING Performs signal processing for MiRS
%   Takes scenario struct as input, retuns scenario.cube struct containing
%   processed Range-Doppler cube

%% Unpack Variables

radarsetup = scenario.radarsetup;
simsetup = scenario.simsetup;

%% Define Constants

c = physconst('LightSpeed');
lambda = c/radarsetup.f_c;

%% Perform Range FFT

% Calculate FFT Size
N_r = 2^ceil(log2(size(scenario.rx_sig,1)));

% Apply windowing
expression = '(size(scenario.rx_sig,1)).*scenario.rx_sig;';
expression = [ radarsetup.r_win, expression];
cube.range_cube = eval(expression);

% FFT across fast time dimension
cube.range_cube = fft(scenario.rx_sig, N_r, 1);

% Remove negative complex frequencies
cube.range_cube = cube.range_cube(1:ceil(end/2),:,:);

%% Perform Doppler FFT

% Calculate FFT size
N_d = 2^ceil(log2(size(cube.range_cube,2)));

% Apply windowing
expression = '(size(cube.range_cube,2))).*cube.range_cube;';
expression = ['transpose(', radarsetup.d_win, expression];
cube.rd_cube = eval(expression);

% Clear range cube
if simsetup.clear_cube
    cube.range_cube = [];
end

% FFT across slow time dimension
cube.rd_cube = fftshift(fft(cube.rd_cube, N_d, 2), 2);

% Wrap max negative frequency and positive frequency
cube.rd_cube(:,(end+1),:) = cube.rd_cube(:,1,:);

%% Perform Angle FFTs

% Calculate FFT size 
if isempty(radarsetup.n_az)
    N_az = 2^ceil(log2(size(cube.rd_cube, 3)));
else
    N_az = radarsetup.n_az;
end


% Apply elevation windowing
if strcmp(radarsetup.az_win, 'none')
    cube.angle_cube = cube.rd_cube;
else
    expression = '(size(cube.rd_cube,3)), [2 3 1]).*cube.rd_cube;';
    expression = ['permute(', radarsetup.az_win, expression];
    cube.angle_cube = eval(expression);
end

% Clear MIMO cube
if simsetup.clear_cube
    cube.rd_cube = [];
end

% FFT across angle dimensions
cube.angle_cube = fftshift(fft(cube.angle_cube, N_az, 3), 3);

% Wrap max negative frequency and positive frequency
cube.angle_cube(:,:,(end+1)) = cube.angle_cube(:,:,1);

%% Calculate Power Cube

% Take square magnitude of radar cube
cube.pow_cube = abs(cube.angle_cube).^2;

% Clear angle cube
if simsetup.clear_cube
    cube.angle_cube = [];
end

%% Derive Axes

% Derive Range axis
cube.range_res = ((size(scenario.rx_sig,1) + radarsetup.drop_s)/N_r)*(c/(2*radarsetup.bw));
cube.range_axis = ((1:(N_r/2))-1)*cube.range_res;

% Derive Doppler axis
cube.vel_res = lambda/(2*radarsetup.t_ch*radarsetup.n_p*radarsetup.n_tx);
cube.vel_axis = ((-N_d/2):(N_d/2))*cube.vel_res;

% Derive Azimuth axis
cube.azimuth_axis = -asind(((-N_az/2):(N_az/2))*(2/N_az));
cube.azimuth_res = min(abs(diff(cube.azimuth_axis)));

end