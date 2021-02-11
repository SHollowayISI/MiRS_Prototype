%% MiRS - Target Initialization File
%{

    Sean Holloway
    MiRS Target Init File
    
    This file specifies the list of targets used in the MiRS simulation
    system.

    
    
%}

%% Target Setup

% Target positions in meters
tgt_pos = [500; ...
            0; ...
            0];
if test_mode
    tgt_pos = tgt_pos_in;
end

% Target velocities in meters per second
tgt_vel = [ 0; ...
            0; ...
            0];     
if test_mode
    tgt_vel = tgt_vel_in;
end

%% Set Target RCS

% Target RCS in dBm^2
tgt_rcs_dbmm = [23];

% Convert to absolute
tgt_rcs = db2pow(tgt_rcs_dbmm);

%% Save Target List Structure

% Target list
scenario.target_list = struct( ...
    'pos',              tgt_pos, ...
    'vel',              tgt_vel, ...
    'rcs',              tgt_rcs);



