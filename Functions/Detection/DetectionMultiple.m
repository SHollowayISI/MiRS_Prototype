function [detection] = DetectionMultiple(scenario)
%DETECTIONMULTIPLE Performs multiple-CPI detection for MiRS project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing information about detected targets.

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
cube = scenario.cube;
sim = scenario.sim;

%% Perform Binary Integration

% Find over-threshold detections
bw_cube = (detection.detect_cube_multi >= radarsetup.det_m);

% Average power for multiple-detection indices
avg_cube = bw_cube .* (detection.pow_cube_multi ./ detection.detect_cube_multi);
avg_cube(isnan(avg_cube)) = 0;

%% Determine Individual Object Coordinates

% Find connected objects in R-D cube
cc = bwconncomp(bw_cube);
regions = regionprops(cc, avg_cube, 'WeightedCentroid');

% Generate list of detection coordinates
detection.detect_list.range = [];
detection.detect_list.vel = [];
detection.detect_list.aoa = [];
detection.detect_list.cart = [];
detection.detect_list.SNR = [];
detection.detect_list.num_detect = length(regions);

% Determine Centroid of azimuth-elevation slice
for n = 1:length(regions)
    
    % Shorten variable for ease of typing
    ind = regions(n).WeightedCentroid;
    
    % Store direct coordinates
    detection.detect_list.range(end+1) = interp1(cube.range_axis, ind(2));
    detection.detect_list.vel(end+1) = interp1(cube.vel_axis, ind(1));
    
    % Store SNR
    detection.detect_list.SNR(end+1) = 10*log10(max(avg_cube(cc.PixelIdxList{n}), [], 'all')) ...
        - detection.noise_pow;
    
    % Find angle of attack using AoA estimator
    ant_slice = squeeze(cube.rd_cube(round(ind(2)), round(ind(1)), :))';
    ang_list = zeros(radarsetup.n_tx,1);
    for tx_num = 1:radarsetup.n_tx
        ind = (tx_num-1)*radarsetup.n_rx + (1:radarsetup.n_rx);
        [~, ang_list(tx_num)] = sim.AoAReal(ant_slice(ind));
    end
    ang = -1 * mean(ang_list);
%     [~, ang] = sim.AoA(ant_slice_corrected);
    
    % Correct velocity-bearing coupling due to TDM-MIMO
    detection.detect_list.aoa(end+1) = ang - radarsetup.v_az_coeff * detection.detect_list.vel(end);
    
    % Store derived coordinates
    detection.detect_list.cart(:,end+1) = detection.detect_list.range(end) * ...
        [cosd(detection.detect_list.aoa(end)); sind(detection.detect_list.aoa(end))];
    
    % Search through previous detections for targets within group radius
    if radarsetup.rm_group
        
        for m = 1:(length(detection.detect_list.range)-1)
            
            if m > (length(detection.detect_list.range)-1)
                break;
            end
            
            % Check distance in each dimension
            dist = abs([detection.detect_list.range(m) - detection.detect_list.range(end), ...
                detection.detect_list.vel(m) - detection.detect_list.vel(end), ...
                detection.detect_list.aoa(m) - detection.detect_list.aoa(end)]);
            
            % If in radius, keep larger SNR of two
            if all(dist < radarsetup.rm_rad)
                if detection.detect_list.SNR(m) < detection.detect_list.SNR(end)
                    detection.detect_list.range(m) = [];
                    detection.detect_list.vel(m) = [];
                    detection.detect_list.aoa(m) = [];
                    detection.detect_list.cart(:,m) = [];
                    detection.detect_list.SNR(m) = [];
                    detection.detect_list.num_detect = detection.detect_list.num_detect - 1;
                else
                    detection.detect_list.range(end) = [];
                    detection.detect_list.vel(end) = [];
                    detection.detect_list.aoa(end) = [];
                    detection.detect_list.cart(:,end) = [];
                    detection.detect_list.SNR(end) = [];
                    detection.detect_list.num_detect = detection.detect_list.num_detect - 1;
                    break;
                end
            end
            
        end
    end
    
end

end

