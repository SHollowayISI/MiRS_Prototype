
% Iterate over multiple values of error_rms to find breakpoint and log cost
% change threshold

var_list = 0:0.001:1;
iterations_per_error = 100;

break_list = zeros(length(var_list), iterations_per_error);
thresh_list = zeros(length(var_list), iterations_per_error);
error_list = zeros(length(var_list), iterations_per_error);

for m = 1:length(var_list)
    
    str = ['List number ', num2str(m)];
    disp(str);
    
    for n = 1:iterations_per_error
        
        error_var = var_list(m);
        MiRS_Multilateration_3Dimension;
        
        break_list(m, n) = breakpoint;
        thresh_list(m, n) = thresh;
        error_list(m, n) = rms(e);
        
    end
end

save('Results.mat', 'break_list', 'thresh_list', 'error_list');