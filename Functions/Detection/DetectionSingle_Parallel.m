function [detection] = DetectionSingle_Parallel(scenario)
%DETECTIONSINGLE_PARALLEL Performs target detection for MiRS project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing information about detected targets. Uses MATLAB's
%   Parallel Computing Toolbox for speedup of CFAR operation.

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
cube = scenario.cube;

%% Perform Detection

% Estimate noise power
noise_pow = pow2db(median(mean(cube.pow_cube, 1), 'all'));
detection.noise_pow = noise_pow;

% Separate out CFAR object for parallel processing
CFAR_Object = scenario.sim.CFAR;

% Set up index map
[n_rng, n_dop, ~] = size(cube.pow_cube);
pad = radarsetup.num_guard + radarsetup.num_train;
rng_ax = (pad(1) + 1):(n_rng-pad(1));
dop_ax = (pad(2) + 1):(n_dop-pad(2));
idx = [];
idx(1,:) = repmat(rng_ax, 1, length(dop_ax));
idx(2,:) = reshape(repmat(dop_ax, length(rng_ax), 1), 1, []);
    
% Separate objects for parallel processing
pow_cube = cube.pow_cube;
detect_cube = zeros(size(pow_cube(:,:,1:(size(pow_cube, 3)-1))));

% Loop across angle slices
parfor angle_slice = 1:(size(cube.pow_cube, 3)-1)
    
    %% Perform CFAR Detection
    
    % Perform CFAR detection
    cfar_out = CFAR_Object(pow_cube(:,:,angle_slice), idx);
    
    % Reshape to radar cube size
    cfar_out = reshape(cfar_out, length(rng_ax), length(dop_ax));
    cfar_out = [zeros(pad(1), length(dop_ax)); cfar_out; zeros(pad(1), length(dop_ax))];
    cfar_out = [zeros(n_rng, pad(2)), cfar_out, zeros(n_rng, pad(2))];
    
    % Perform image dilation
    %             if radarsetup.dilate
    %                 str = ['se = strel(radarsetup.dilate_shape, ', radarsetup.dilate_args, ');'];
    %                 eval(str);
    %                 cfar_out = imdilate(cfar_out, se);
    %             end
    
    % Save detection cube
    detect_cube(:,:,angle_slice) = cfar_out;
    
end

% Clean parallel ends
detection.detect_cube = detect_cube;

% Wrap ends of angle FFT
detection.detect_cube(:,:,size(cube.pow_cube, 3)) = detection.detect_cube(:,:,1);

%% Update Multiple CPI List

% Initialize multi-frame arrays if not created
if isempty(detection.detect_cube_multi)
    detection.detect_cube_multi = zeros(size(detection.detect_cube));
    detection.pow_cube_multi = zeros(size(cube.pow_cube));
end

% Add to number of detections per cell
detection.detect_cube_multi = detection.detect_cube_multi + detection.detect_cube;

% Add to power cube
detection.pow_cube_multi = detection.pow_cube_multi + cube.pow_cube .* detection.detect_cube;

end



