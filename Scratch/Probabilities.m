%% Setup

SNR_min = 9.7615;
P_d_req = 0.95;
k = 1;
% P_k = 0.06;
num_samp = (45e6 * 1.6e-3) - 7200;
n_r = 2^ceil(log2(num_samp));
N = n_r * 64 * 8;

%% Workspace

% Pfa from P_k
% [P_d_mn, P_d, P_fa, P_fa_mn] = Pd_from_Pk(P_k, k, SNR_min, N, 2, 5);

% P_k from Pd
P_k = 0.01:0.01:0.25;
for b = 1:length(P_k)
    [Pd_mn_out(b), Pd_out(b), Pfa_out(b), Pfa_mn_out(b)] = Pd_from_Pk(P_k(b), k, SNR_min, N, 2, 5);
end

[P_k', Pd_mn_out', Pd_out', Pfa_out', Pfa_mn_out']

%% User defined functions

function [z] = invprob(p)
    z = sqrt(2) * erfinv(2*p - 1);
end

function P_fa_mn = Pfa_from_K(k, n, z)
    
%     P_fa_plus = ((2*k*n + z*z*n) + sqrt((2*k*n + z*z*n)^2 - 4*k*k*(n*n + z*z*n))) / ...
%         (2*(n*n + z*z*n));
    P_fa_mn = ((2*k*n + z*z*n) - sqrt((2*k*n + z*z*n)^2 - 4*k*k*(n*n + z*z*n))) / ...
        (2*(n*n + z*z*n));
    
end

function P_fa_mn = Pfa_Binary(P_fa_single, m, n)
    
    P_fa_mn = zeros(size(P_fa_single));
    for i = m:n
        P_fa_mn = P_fa_mn + nchoosek(n, i) * (P_fa_single).^i .* (1-P_fa_single).^(n-i);
    end

end

function P_fa_single = Pfa_Reverse(P_fa_mn, m, n)

    P_fa_single_in = 10.^(-1:-0.001:-10)';
    P_fa_single_out = Pfa_Binary(P_fa_single_in, m, n);
    
    ind = find(P_fa_single_out < P_fa_mn, 1);
    P_fa_single = sqrt(prod(P_fa_single_in([ind-1, ind])));

end

function P_d_mn = Pd_Binary(P_d_single, m, n)
    
    P_d_mn = zeros(size(P_d_single));
    for i = m:n
        P_d_mn = P_d_mn + nchoosek(n, i) * (P_d_single).^i .* (1-P_d_single).^(n-i);
    end
    
end

function P_d_single = Pd_Reverse(P_d_mn, m, n)

    P_d_single_in = (0:0.001:1)';
    P_d_single_out = Pd_Binary(P_d_single_in, m, n);
    
    ind = find(P_d_single_out > P_d_mn, 1);
    P_d_single = sqrt(prod(P_d_single_in([ind-1, ind])));

end

function P_fa_single = Pfa_from_Pd(P_d_single, SNR)
    
    [d, fa] = rocsnr(SNR);
    close all;
    ind = find(d > P_d_single, 1);
    P_fa_single = sqrt(prod(fa([ind-1, ind])));

end

function P_d_single = Pd_from_Pfa(P_fa_single, SNR)
    
    [d, fa] = rocsnr(SNR);
    close all;
    ind = find(fa > P_fa_single, 1);
    P_d_single = sqrt(prod(d([ind-1, ind])));

end

function [P_d_mn, P_d_single, P_fa_single, P_fa_mn] = Pd_from_Pk(P_k, k, SNR, N, m, n)
    
    P_fa_mn = Pfa_from_K(k, N, invprob(1 - P_k));
    P_fa_single = Pfa_Reverse(P_fa_mn, m, n);
    P_d_single = Pd_from_Pfa(P_fa_single, SNR);
    P_d_mn = Pd_Binary(P_d_single, m, n);

end





