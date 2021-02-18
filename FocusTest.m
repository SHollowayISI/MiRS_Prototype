
%% Bookkeeping

close all;
clear variables;
tic;

%% Setup

% Signal properties
centerFrequencies = [0.5, 1, 3, 6] * 1e9;
bandwidth = 10e6;

% Coordinate setup
numberOfNodes = 10;
heightAboveTerrain = 1e3;
heightVariation = 10;
fieldSize = 1.6e3;
groundDistanceToBase = 10e3;

% Coordinate calculations
baseCoordinates = [groundDistanceToBase; 0; 0];
masterCoordinates = [0; 0; heightAboveTerrain];
slaveCoordinates = [fieldSize * (rand(2, numberOfNodes-1)-0.5);...
    heightAboveTerrain + heightVariation*(rand(1, numberOfNodes-1)-0.5)];
nodeCoordinates = [masterCoordinates, slaveCoordinates];

%% Focusing Algorithm

[tof, td] = CalculateDelay(baseCoordinates, nodeCoordinates)

%% Function declarations

function [timeOfFlight, timeDelay] = CalculateDelay(baseCoord, nodeCoords)

ranges = sqrt(sum((nodeCoords - baseCoord).^2, 1));

timeOfFlight = ranges / physconst('Lightspeed');
maxTOF = max(timeOfFlight);
timeDelay = timeOfFlight - maxTOF;

end


















