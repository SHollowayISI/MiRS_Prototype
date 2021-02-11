
%% To Do:
%{
    - Elimination of branches
%}

%% Bookkeeping

clear variables;
close all;

%% Setup

% Properties
num_nodes = 10;
field_size = 1000;
error_var = 0.01;

% Iteration properties
remove_branches = false;
branch_per_node = 4;
max_iterations = 1000;
cost_diff_thresh = -Inf;

% Calculations
num_var = 3*num_nodes - 6;
num_eq = (num_nodes - 1) * num_nodes / 2;

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = field_size * (rand(num_var + 6, 1)-0.5);
x_guess{1} = field_size * rand(num_var, 1);

% Calculate exact ranges (squared)
n_ind = 3*floor((n-1)) + (1:3) - 6;
m_ind = 3*floor((m-1)) + (1:3) - 6;
r_sq = sum((x_real(n_ind + 6) - x_real(m_ind + 6)).^2, 2);

% Introduce error into range measurement
r_sq = r_sq + error_var*rand(size(r_sq));

%% Branch removal

% Remove random branches from tree
branch_rm = num_nodes - 1 - branch_per_node;
null_ind = zeros(num_nodes, branch_rm);

null_branches = [];
if (remove_branches && (branch_rm > 0))
    for i = 1:num_nodes
        
        others = setdiff((1:num_nodes), i);
        neighbors = mod(i + (-branch_per_node/2:branch_per_node/2) - 1, num_nodes) + 1;
        
        null_ind(i,:) = setdiff(others, neighbors);
        
    end
    
    for i = 1:num_nodes
        for k = 1:branch_rm
            
            null_num = (i - 2) * (i - 1) * 0.5 + null_ind(i, k);
            null_branches = [null_branches, null_num];
            
        end
    end
end


%% Initial step

% Evaluate cost function
s = zeros(length(n),1);
for i = 1:length(n)
    if (n(i) == 2 && m(i) == 1)
        s(i) = x_guess{1}(1)^2 - r_sq(i);
    elseif (n(i) == 3 && m(i) == 1)
        s(i) = sum(x_guess{1}(2:3).^2) - r_sq(i);
    elseif (m(i) == 1)
        s(i) = sum(x_guess{1}(n_ind(i,1:3)).^2) - r_sq(i);
    elseif (n(i) == 3 && m(i) == 2)
        s(i) = (x_guess{1}(2) - x_guess{1}(1))^2 + x_guess{1}(3)^2 - r_sq(i);
    elseif (m(i) == 2)
        s(i) = (x_guess{1}(n_ind(i,1)) - x_guess{1}(1))^2 + sum(x_guess{1}(n_ind(i,2:3)).^2) - r_sq(i); 
    elseif (m(i) == 3)
        s(i) = sum((x_guess{1}(n_ind(i,1:2)) - x_guess{1}(2:3)).^2) + x_guess{1}(n_ind(i,3))^2 - r_sq(i);
    else
        s(i) = sum((x_guess{1}(n_ind(i,:)) - x_guess{1}(m_ind(i,:))).^2) - r_sq(i);
    end
end

% Remove nullified branches
s(null_branches) = [];

cost(1) = sum(s.^2);

%% Main loop

for iter = 2:max_iterations
    
    % Jacobian
    J = zeros(num_eq, num_var);
    for i = 1:num_eq
        if (n(i) == 2 && m(i) == 1)
            J(i, 1) = 2 * (x_guess{iter-1}(1));
        elseif (n(i) == 3 && m(i) == 1)
            J(i, 2:3) = 2 * (x_guess{iter-1}(2:3));
        elseif (m(i) == 1)
            J(i, n_ind(i, :)) = 2 * (x_guess{iter-1}(n_ind(i,:)));
        elseif (n(i) == 3 && m(i) == 2)
            J(i, 2) = 2 * (x_guess{iter-1}(2) - x_guess{iter-1}(1));
            J(i, 3) = 2 * (x_guess{iter-1}(3));
            J(i, 1) = 2 * (x_guess{iter-1}(1) - x_guess{iter-1}(2));
        elseif (m(i) == 2)
            J(i, n_ind(i, 1))   = 2 * (x_guess{iter-1}(n_ind(i,1)) - x_guess{iter-1}(1));
            J(i, n_ind(i, 2:3)) = 2 * (x_guess{iter-1}(n_ind(i,2:3)));
            J(i, 1)             = 2 * (x_guess{iter-1}(1) - x_guess{iter-1}(n_ind(i,2))); 
        elseif (m(i) == 3)
            J(i, n_ind(i, 1:2)) = 2 * (x_guess{iter-1}(n_ind(i,1:2)) - x_guess{iter-1}(2:3));
            J(i, n_ind(i, 3))   = 2 * (x_guess{iter-1}(n_ind(i,3)));
            J(i, 2:3)           = 2 * (x_guess{iter-1}(2:3) - x_guess{iter-1}(n_ind(i,1:2)));
        else
            J(i, n_ind(i,:)) = 2 * (x_guess{iter-1}(n_ind(i,:)) - x_guess{iter-1}(m_ind(i,:)));
            J(i, m_ind(i,:)) = 2 * (x_guess{iter-1}(m_ind(i,:)) - x_guess{iter-1}(n_ind(i,:)));
        end
    end
    
    % Nullify branches
    J(null_branches, :) = [];
    
    % Gauss-Newton step
    x_guess{iter} = x_guess{iter-1} - inv(J' * J) * J' * s;
    
    % Sum of squares
    s = zeros(length(n),1);
    for i = 1:length(n)
        if (n(i) == 2 && m(i) == 1)
            s(i) = x_guess{iter}(1)^2 - r_sq(i);
        elseif (n(i) == 3 && m(i) == 1)
            s(i) = sum(x_guess{iter}(2:3).^2) - r_sq(i);
        elseif (m(i) == 1)
            s(i) = sum(x_guess{iter}(n_ind(i,1:3)).^2) - r_sq(i);
        elseif (n(i) == 3 && m(i) == 2)
            s(i) = (x_guess{iter}(2) - x_guess{iter}(1))^2 + x_guess{iter}(3)^2 - r_sq(i);
        elseif (m(i) == 2)
            s(i) = (x_guess{iter}(n_ind(i,1)) - x_guess{iter}(1))^2 + sum(x_guess{iter}(n_ind(i,2:3)).^2) - r_sq(i);
        elseif (m(i) == 3)
            s(i) = sum((x_guess{iter}(n_ind(i,1:2)) - x_guess{iter}(2:3)).^2) + x_guess{iter}(n_ind(i,3))^2 - r_sq(i);
        else
            s(i) = sum((x_guess{iter}(n_ind(i,:)) - x_guess{iter}(m_ind(i,:))).^2) - r_sq(i);
        end
    end
    
    % Remove nullified branches
    s(null_branches) = [];
    
    cost(iter) = sum(s.^2);
    
    cost_diff = log10(abs(cost(iter) - cost(iter-1)));
    if (cost_diff < cost_diff_thresh)
        iter
        break;
    end
    
end

%% Error estimation & Plotting

% Degree of freedom reduction and error estimation
x_result = reshape([0; 0; 0; x_guess{end}(1); 0; 0; x_guess{end}(2:3); 0; x_guess{end}(4:end)], ...
    3, num_nodes)';
x_exact = reshape(x_real, 3, num_nodes)';
x_exact = x_exact - x_exact(1,:);

% Perform rotation
%
x1 = x_exact(2,1);
y1 = x_exact(2,2);
x2 = x_result(2,1);

% Rotate p2 into x-z plane (around z)
theta = atan2d(-y1*x2, x1*x2);
R = [cosd(theta), -sind(theta), 0; sind(theta), cosd(theta), 0; 0, 0, 1];
x_rotated = transpose(R * x_exact');

% Rotate p2 onto x-axis (around y)
x1 = x_rotated(2,1);
y1 = x_rotated(2,3);
x2 = x_result(2,1);

theta = -atan2d(-y1*x2, x1*x2);
R = [cosd(theta), 0, sind(theta); 0, 1, 0; -sind(theta), 0, cosd(theta)];
x_rotated = transpose(R * x_rotated');

% Rotate p3 into place (around x)
x1 = x_rotated(3,2);
y1 = x_rotated(3,3);
x2 = x_result(3,2);

theta = atan2d(-y1*x2, x1*x2);
R = [1, 0, 0; 0, cosd(theta), -sind(theta); 0, sind(theta), cosd(theta)];
x_rotated = transpose(R * x_rotated');

% Perform reflection
flip = sign(x_rotated(4,3))/sign(x_result(4,3));
x_rotated = x_rotated .* [1, 1, flip];

% Calculate error
e = sum((x_rotated - x_result).^2, 2);
rms(e)

%% Plotting
close all;

% Plot points
%
figure;
scatter3(x_result(:,1), x_result(:,2), x_result(:,3), '.', 'r');
% hold on;
% scatter3(x_result(2,1), x_result(2,2), x_result(2,3), '.', 'b');
% hold on;
% scatter3(x_result(3,1), x_result(3,2), x_result(3,3), '.', 'g');
hold on;
scatter3(x_rotated(:,1), x_rotated(:,2), x_rotated(:,3));
% hold on;
% scatter3(x_rotated(2,1), x_rotated(2,2), x_rotated(2,3));
% hold on;
% scatter3(x_rotated(3,1), x_rotated(3,2), x_rotated(3,3));
grid on;
xlim([-1 1] * field_size)
ylim([-1 1] * field_size)
zlim([-1 1] * field_size)
%}

% Cost function over iterations
log_cost_delta = log10(abs(diff(cost)));
% thresh = max(log_cost_delta((end-100):end))

% breakpoint = find(log_cost_delta < thresh, 1)

%{
figure;
plot(log_cost_delta);
% hold on;
% plot(ones(size(cost))*thresh);
grid on;
%}







