%% MiRS - Main Processing Loop
%{

    Sean Holloway
    MiRS Main Processing Loop
    
    This file specifies performs simulation, signal processing, detection,
    data processing, and results collection for MiRS system.
    
%}

%% Main Loop

% Start timing for estimation
timeStart(scenario);

for loop = 1:scenario.simsetup.num_frames
    
    scenario.flags.frame = loop;
    
    for cpi = 1:scenario.radarsetup.cpi_fr
        
        scenario.flags.cpi = cpi;
        
        %% Radar Simulation
        
        % Run simulation to retrieve fast time x slow time x rx-channel signal
        if scenario.simsetup.par_sim
            scenario = RadarSimulation_Parallel_StaticTarget(scenario);
%             scenario = RadarSimulation_Parallel(scenario);
        else
            scenario = RadarSimulation(scenario);
        end
        
        %% Signal Processing
        
        % Perform signal processing on received signal
        scenario.cube = SignalProcessing(scenario);
        
        %% Single CPI Data Processing
        
        % Perform single-frame radar detection
        if scenario.simsetup.par_cfar
            scenario.detection = DetectionSingle_Parallel(scenario);
        else
            scenario.detection = DetectionSingle(scenario);
        end
        
        %% Loop Update Procedures
        
        % Read out CPI update
        CPIUpdate(scenario);
        
        % Read out estimated time of completion
        timeUpdate(scenario, 1, 'loops')
        
    end
    
    %% Multiple CPI Data Processing
    
    % Perform binary integration and coordinate determination
    scenario.detection = DetectionMultiple(scenario);
    
    % Read out detection data
    if scenario.simsetup.readout
        readOut(scenario);
    end
    
    % Save detection data
    saveMulti(scenario);
    
    % Update multi-target tracking system
%     scenario.multi = Tracking(scenario);
    
end

%% Multi Frame Data Visualization

% Display result visualization plots
% viewTracking(scenario);








