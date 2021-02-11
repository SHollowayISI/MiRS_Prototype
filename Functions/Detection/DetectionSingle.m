function [detection] = DetectionSingle(scenario)
%DETECTIONSINGLE Performs target detection for MiRS project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing information about detected targets.

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
cube = scenario.cube;
flags = scenario.flags;

%% Perform Detection

% Estimate noise power
detection.noise_pow = pow2db(median(mean(cube.pow_cube, 1), 'all'));

% Loop across angle slices
for angle_slice = 1:(size(cube.pow_cube, 3)-1)
    
    rd_cube = cube.pow_cube(:,:,angle_slice);
    
    switch radarsetup.detect_type
        case 'threshold'
            %% Perform Threshold Detection
            
            % Calculate threshold in absolute
            abs_thresh = db2pow(radarsetup.thresh + detection.noise_pow);
            
            % Perform detection
            detection.detect_cube(:,:,angle_slice) = (rd_cube > abs_thresh);
            
        case 'CFAR'
            %% Perform CFAR Detection
            
            % Set up index map
            [n_rng, n_dop] = size(rd_cube);
            pad = radarsetup.num_guard + radarsetup.num_train;
            rng_ax = (pad(1) + 1):(n_rng-pad(1));
            dop_ax = (pad(2) + 1):(n_dop-pad(2));
            
            idx = [];
            idx(1,:) = repmat(rng_ax, 1, length(dop_ax));
            idx(2,:) = reshape(repmat(dop_ax, length(rng_ax), 1), 1, []);
            
            % Perform CFAR detection
            cfar_out = scenario.sim.CFAR(rd_cube, idx);
            
            % Reshape to radar cube size
            cfar_out = reshape(cfar_out, length(rng_ax), length(dop_ax));
            cfar_out = [zeros(pad(1), length(dop_ax)); cfar_out; zeros(pad(1), length(dop_ax))];
            cfar_out = [zeros(n_rng, pad(2)), cfar_out, zeros(n_rng, pad(2))];
            
            % Perform image dilation
            if radarsetup.dilate
                str = ['se = strel(radarsetup.dilate_shape, ', radarsetup.dilate_args, ');'];
                eval(str);
                cfar_out = imdilate(cfar_out, se);
            end
            
            % Save detection cube
            detection.detect_cube(:,:,angle_slice) = cfar_out;
            
    end
end

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



