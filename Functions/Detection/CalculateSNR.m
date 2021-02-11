function [SNR] = CalculateSNR(scenario,RCS,Range)
%CALCULATESNR Calculates SNR of target for MiRS project
%   Takes radar scenario object, RCS (absolute) value, and Range 
%   as input, provides SNR value as output.

%% Unpack Variables
radarsetup = scenario.radarsetup;

%% Calculate SNR

n_tx = radarsetup.n_tx;
n_rx = radarsetup.n_rx;

lambda = physconst('LightSpeed')/radarsetup.f_c;
t_cpi = radarsetup.t_ch*radarsetup.n_p*n_tx;
Lr = db2pow(radarsetup.rf_sys_loss);
NF = db2pow(radarsetup.rx_nf);
c = (4*pi)^3;
n = physconst('Boltzmann')*290;

Gt = radarsetup.ant_gain;
Gr = Gt;
total_pow = radarsetup.tx_pow;

SNR_abs = (total_pow * Gt * Gr * t_cpi * lambda * lambda * RCS) ...
    ./ (c * (Range.^4) * n * NF * Lr);

SNR = pow2db(SNR_abs);

end

