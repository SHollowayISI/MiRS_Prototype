% ClassDef File for MiRS Radar Scenario

classdef RadarScenario < handle
    properties
        target_list
        simsetup
        radarsetup
        sim
        rx_sig
        cube
        detection
        flags
        timing
        results
        multi
    end
    
    methods
        
        function RadarScenario = RadarScenario()
            % Initialize structure of target list
            RadarScenario.target_list = struct( ...
                'pos',      [], ...
                'vel',      [], ...
                'rcs',      []);
            
            RadarScenario.multi.detect_list = {};
            
            RadarScenario.multi.track_list = {};
            
            RadarScenario.multi.active_tracks = [];
            
            RadarScenario.detection.detect_cube_multi = [];
        end
        
        function timeStart(RadarScenario)
            % Begin timing for progress readout
            RadarScenario.timing.timing_logical = true;
            RadarScenario.timing.startTime = tic;
            RadarScenario.timing.TimeDate = now;
            RadarScenario.timing.numLoops = ...
                RadarScenario.simsetup.num_frames * ...
                RadarScenario.radarsetup.cpi_fr;
            RadarScenario.timing.timeGate = 0;
        end
        
        function timeUpdate(RadarScenario, repetition, rep_method)
            
            if ~RadarScenario.timing.timing_logical
                error('Must use method timeStart() before timeUpdate()');
            end
            
            % Calculate progress through simulation
            loops_complete = (RadarScenario.flags.frame-1)*RadarScenario.radarsetup.cpi_fr + ...
                RadarScenario.flags.cpi;
            percent_complete = 100*loops_complete/RadarScenario.timing.numLoops;
            
            % Calculate remaining time in simulation
            nowTimeDate = now;
            elapsedTime = nowTimeDate - RadarScenario.timing.TimeDate;
            estComplete = nowTimeDate + ((100/percent_complete)-1)*elapsedTime;
            
            % Form message to display in command window
            if loops_complete > 1
                message_l = sprintf('%d Loops complete out of %d', loops_complete, RadarScenario.timing.numLoops);
            else
                message_l = sprintf('%d Loop complete out of %d', loops_complete, RadarScenario.timing.numLoops);
            end
            message_p = [sprintf('Percent complete: %0.0f', percent_complete), '%'];
            message_t = ['Estimated time of completion: ', datestr(estComplete)];
            
            % Display current progress
            switch rep_method
                case 'loops'
                    if (mod(loops_complete, repetition) == 1) || (repetition == 1)
                        disp(message_l);
                        disp(message_p);
                        disp(message_t);
                        disp('');
                    end
                    
                case 'time'
                    if ((RadarScenario.timing.timeGate == 0) || ...
                            (toc > repetition + RadarScenario.timing.timeGate))
                        disp(message_p);
                        disp(message_t);
                        disp('');
                        RadarScenario.timing.timeGate = toc;
                    end
                    
            end
        end
        
        function CPIUpdate(RadarScenario)
            fprintf('\nCPI %d complete out of %d per frame.\n', ...
                RadarScenario.flags.cpi, ...
                RadarScenario.radarsetup.cpi_fr);
        end
        
        function readOut(RadarScenario)

            num_detect = RadarScenario.detection.detect_list.num_detect;
            
            if num_detect > 0
                if num_detect > 1
                    fprintf('\n%d Targets Detected:\n', num_detect);
                else
                    fprintf('\n%d Target Detected:\n', num_detect);
                end
                
                for n = 1:num_detect
                    fprintf('\nTarget #%d Coordinates:\n', n);
                    fprintf('Range: %0.2f [m]\n', ...
                        RadarScenario.detection.detect_list.range(n));
                    fprintf('Velocity: %0.2f [m/s]\n', ...
                        RadarScenario.detection.detect_list.vel(n));
                    fprintf('Bearing: %0.2f [deg]\n', ...
                        RadarScenario.detection.detect_list.aoa(n));
                    fprintf('SNR: %0.1f [dB]\n', ...
                        RadarScenario.detection.detect_list.SNR(n));
                end
            else
                fprintf('\nNo Targets Detected\n');
            end
            
            if length(RadarScenario.target_list.rcs) == 1
                idealSNR = CalculateSNR(RadarScenario, RadarScenario.target_list.rcs, ...
                    sqrt(sum(RadarScenario.target_list.pos.^2)));
                fprintf('\nIdeal SNR of Target: %0.1f\n', idealSNR);
            end
            
        end
        
        function saveMulti(RadarScenario)
            
            RadarScenario.multi.detect_list{RadarScenario.flags.frame} = ...
                RadarScenario.detection.detect_list;
            
        end
        
        function viewTargets(RadarScenario)
            % Find axis limits
            max_dist = sqrt(max(sum(RadarScenario.target_list.pos(1:2,:).^2, 1)));
            % Show 3-D scatter plot of target locations and transceiver
            figure('Name', 'Target Scatter Plot')
            scatter(RadarScenario.target_list.pos(2,:), ...
                RadarScenario.target_list.pos(1,:))
            hold on;
            scatter(RadarScenario.simsetup.radar_pos(2), ...
                RadarScenario.simsetup.radar_pos(1), '.', 'r');
            grid on;
            xlim(max_dist * [-1 1])
            ylim(max_dist * [0 2])
            title('Target Scatter Plot')
            xlabel('Cross-Range Distance [m]','FontWeight','bold')
            ylabel('Down-Range Distance [m]','FontWeight','bold')
        end
        
        function viewRDCube(RadarScenario, angleSlice, graphType)
            if strcmp(graphType, 'heatmap')
                figure('Name', 'Range-Doppler Heat Map');
                imagesc(RadarScenario.cube.vel_axis, ...
                    RadarScenario.cube.range_axis, ...
                    10*log10(RadarScenario.cube.pow_cube(:,:,angleSlice)))
                title('Range-Doppler Heat Map')
                set(gca,'YDir','normal')
                xlabel('Velocity [m/s]','FontWeight','bold')
                ylabel('Range [m]','FontWeight','bold')
            else
                figure('Name', 'Range-Doppler Surface');
                surf(RadarScenario.cube.vel_axis, ...
                    RadarScenario.cube.range_axis, ...
                    10*log10(RadarScenario.cube.pow_cube(:,:,angleSlice)), ...
                    'EdgeColor', 'none')
                title('Range-Doppler Surface')
                set(gca,'YDir','normal')
                xlabel('Velocity [m/s]','FontWeight','bold')
                ylabel('Range [m]','FontWeight','bold')
                zlabel('FFT Log Intensity [dB]','FontWeight','bold')
            end
            
        end
        
        function viewRACube(RadarScenario, dopplerSlice, graphType)
            switch graphType
                case 'heatmap'
                    figure('Name', 'Range-Azimuth Heat Map');
                    imagesc(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        10*log10(squeeze(RadarScenario.cube.pow_cube(:,dopplerSlice,:))))
                    title('Range-Azimuth Heat Map')
                    set(gca,'YDir','normal')
                    xlabel('Azimuth Angle [degree]','FontWeight','bold')
                    ylabel('Range [m]','FontWeight','bold')
                case 'surface'
                    figure('Name', 'Range-Azimuth Surface');
                    surf(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        10*log10(squeeze(RadarScenario.cube.pow_cube(:,dopplerSlice,:))), ...
                        'EdgeColor', 'none')
                    title('Range-Azimuth Surface')
                    xlabel('Azimuth Angle [degree]','FontWeight','bold')
                    ylabel('Range [m]','FontWeight','bold')
                    zlabel('FFT Log Intensity [dB]','FontWeight','bold')
                case 'PPI'
                    % Generate coordinate grid
                    
                    % Plot PPI along coordinate grid
                    figure('Name', 'Range-Azimuth PPI');
                    surf(RadarScenario.cube.x_grid, ...
                        RadarScenario.cube.y_grid, ...
                        10*log10(squeeze(RadarScenario.cube.pow_cube(:,dopplerSlice,:))), ...
                        'EdgeColor', 'none')
                    title('Range-Azimuth PPI')
                    xlabel('Cross-Range Distance [m]','FontWeight','bold')
                    ylabel('Down-Range Distance [m]','FontWeight','bold')
                    zlabel('FFT Log Intensity [dB]','FontWeight','bold')
            end
            
        end
        
        function viewDetections(RadarScenario, graphType)
            switch graphType
                case 'heatmap'
                    figure('Name', 'Detection Heatmap')
                    imagesc(RadarScenario.cube.azimuth_axis, ...
                        RadarScenario.cube.range_axis, ...
                        squeeze(any((RadarScenario.detection.detect_cube_multi >= RadarScenario.radarsetup.det_m), 2)))
                    set(gca, 'YDir', 'Normal')
                    title('Detection Heatmap')
                    xlabel('Bearing [degree]', 'FontWeight', 'bold')
                    ylabel('Ragne [m]', 'FontWeight', 'bold')
                case 'PPI'
                    % Generate coordinate grid
                    x_grid = RadarScenario.cube.range_axis' .* cosd(RadarScenario.cube.azimuth_axis);
                    y_grid = RadarScenario.cube.range_axis' .* sind(RadarScenario.cube.azimuth_axis);
                    % Plot PPI along coordinate grid
                    figure('Name', 'Detection PPI');
                    surf(x_grid, y_grid, ...
                        int8(squeeze(any((RadarScenario.detection.detect_cube_multi >= RadarScenario.radarsetup.det_m), 2))), ...
                        'EdgeColor', 'none')
                    view(270,90)
                    title('Detection PPI')
                    xlabel('Cross-Range Distance [m]','FontWeight','bold')
                    ylabel('Down-Range Distance [m]','FontWeight','bold')
                case 'scatter'
                    % Find axis limits
                    max_dist = sqrt(max(sum(RadarScenario.detection.detect_list.cart.^2, 1)));
                    % Plot scatter plot of detections
                    figure('Name', 'Detection Scatter Plot')
                    scatter(RadarScenario.detection.detect_list.cart(2,:), ...
                        RadarScenario.detection.detect_list.cart(1,:));
                    hold on;
                    scatter(RadarScenario.simsetup.radar_pos(2), ...
                        RadarScenario.simsetup.radar_pos(1), '.', 'r');
                    grid on;
                    title('Detection Scatter Plot')
                    xlim(max_dist * [-1 1])
                    ylim(max_dist * [0 2])
                    xlabel('Cross-Range Distance [m]','FontWeight','bold')
                    ylabel('Down-Range Distance [m]','FontWeight','bold')
            end
        end
        
        function viewTracking(RadarScenario)
            % Pass in variables
            track_list  = RadarScenario.multi.track_list;
            % Generate plot
            figure('Name', 'Tracking Results Scatter Plot');
            % Add tracks to plot
            for n = 1:length(track_list)
                % Scatter plot if false alarm
                if track_list{n}.false_alarm
                    scatter3(track_list{n}.det_list(1,:), ...
                        track_list{n}.det_list(2,:), ...
                        track_list{n}.det_list(3,:), ...
                        'r', '+');
                    hold on;
                % Line of track if not
                else
                    scatter3(track_list{n}.det_list(1,:), ...
                        track_list{n}.det_list(2,:), ...
                        track_list{n}.det_list(3,:), ...
                        'k', '+');
                    hold on;
                    plot3(track_list{n}.est_list(1,:), ...
                        track_list{n}.est_list(3,:), ...
                        track_list{n}.est_list(5,:), ...
                        'g');
                    hold on;
                end
            end
            % Add radar location to plot
            scatter3(0, 0, 0, 'filled', 'r');
            % Correct plot limits
            ax = gca;
            ax.YLim = [-ax.XLim(2)/2, ax.XLim(2)/2];
            ax.ZLim = [-ax.XLim(2)/2, ax.XLim(2)/2];
            % Add labels
            xlabel('Down Range Distance [m]', 'FontWeight', 'bold')
            ylabel('Cross Range Distance [m]', 'FontWeight', 'bold')
            zlabel('Altitude [m]', 'FontWeight', 'bold')
        end
             
    end
end






