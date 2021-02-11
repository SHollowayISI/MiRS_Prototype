
clear variables;
close all;
tic;

ranges = 50:50:1000;
iterations = 10;
snr_list = zeros(length(ranges), iterations);
calc_list = snr_list;

for n = 1:length(ranges)
    
    for m = 1:iterations
        
        range_in = ranges(n);
        tgt_pos = range_in * [1; 0; 0];
        
        FullSystem
        
        calc_list(n, m) = idealSNR;
        if(scenario.multi.detect_list{1}.num_detect > 0)
            [~, ind] = min(abs(scenario.multi.detect_list{1}.range - range_in));
            snr_list(n, m) = scenario.multi.detect_list{1}.SNR(ind);
        else
            snr_list(n, m) = nan;
        end
        scenario = [];
        
    end
    
end

plot(snr_list)
hold on;
plot(calc_list)