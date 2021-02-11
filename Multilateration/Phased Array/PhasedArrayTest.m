lambda = 1;
array_center = [0; 0; 100];
array_diff = [0.5; 0; 0];
num_elements = 10;

array = array_center + (array_diff .* ((1:num_elements) - mean(1:num_elements)));

ground_pt = [0; 0; 0];
ranges = rangeangle(array, ground_pt);

steering_angle = 0;
phase_offset = exp(-1i * (1:num_elements) * pi * sind(steering_angle))';

phases = exp(-1i * 2 * pi * ranges / lambda)';
phases = phases ./ phases(ceil(end/2));
phases = phases .* phase_offset;

range_factor = (ranges') .^ -4;
range_factor = range_factor / median(range_factor);

pattern_factor = 

power = 20*log10(abs(sum(phases .* range_factor)))