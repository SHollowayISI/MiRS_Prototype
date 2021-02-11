function [scenario_out] = RadarSimulation(scenario)
%RADARSIMULATION Generates simulated radar response for MiRS
%   Takes scenario.sim, .simsetup, .radarsetup, .target_list as inputs and
%   outputs scenario.rx_sig struct containing the received signal.

%% Unpack Variables

scenario_out = scenario;
target_list = scenario.target_list;
sim = scenario.sim;
simsetup = scenario.simsetup;
radarsetup = scenario.radarsetup;
flags = scenario.flags;

%% Setup

% Physical constants
c = physconst('LightSpeed');

% Derived variables
num_ant = radarsetup.n_rx * radarsetup.n_tx;

%% Simulation

% Allocate size for signals
rx_sig = zeros(radarsetup.n_s - radarsetup.drop_s, radarsetup.n_p, num_ant);

% Generate the pulse
base_sig = sim.waveform();

% Transmit the pulse
tx_sig = sim.transmitter(base_sig);

% NOTE: THE FOLLOWING TWO SETTINGS HAVE BEEN REMOVED FROM PARFOR FOR STATIC TARGET
% Update target position and velocities
[target_list.pos, target_list.vel] = sim.target_plat(radarsetup.t_fr);

% Get range and angle to all targets
[~, tgt_ang] = rangeangle(target_list.pos, simsetup.radar_pos);

% Set up signals to the point where stochastic noise is added
for chirp = 1:radarsetup.n_tx
    
    % Set up weights to implement TDM-MIMO
    weights{chirp} = [0; 0];
    weights{chirp}(chirp) = 1;
    
    % Radiate signal towards all targets
    sig = sim.radiator(tx_sig, tgt_ang, weights{chirp});
    
    % Propogate signal to the target through two-way channel
    sig = sim.channel(sig, ...
        simsetup.radar_pos, target_list.pos, ...
        simsetup.radar_vel, target_list.vel);
    
    % Dechirp signal (Done early to reduce processing overhead)
    sig = dechirp(sig, base_sig);
    
    % Decimate signal to baseband
    decim_sig{chirp} = decimate(sig, radarsetup.sim_rate);
    
    % Reflect the pulse off of the target
    decim_sig{chirp} = sim.target(decim_sig{chirp});
    
    % Collect the reflected target at the antenna array
    decim_sig{chirp} = sim.collector(decim_sig{chirp}, tgt_ang);
    
end


parfor block = 1:radarsetup.n_p
    
    for chirp = 1:radarsetup.n_tx
        
        % Apply receiver noise and gains
        par_sig = sim.receiver(decim_sig{chirp});
        
        % Apply phase coding to result signal
        par_sig = reshape(par_sig, length(par_sig), 1, []);
        
        % Calculate indices for virtual antenna
        ch_ind = (1:radarsetup.n_rx) + radarsetup.n_rx * (chirp - 1);
        
        % Save Rx signal by fast time x slow time x virtual antenna
        par_rx_sig{block}(:,1,ch_ind) = par_sig((radarsetup.drop_s+1):end,:,:);
        
    end
    
end

for block = 1:radarsetup.n_p
    
    rx_sig(:,block,:) = par_rx_sig{block};
    
end


%% Re-pack Variables

scenario_out.rx_sig = rx_sig;


end