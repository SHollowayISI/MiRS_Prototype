
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
branch_per_node = 6;
max_iterations = 1000;
cost_diff_thresh = -10;
cost_thresh = 1;
lambda_init = 0.01;
nu = 1.1;

% Calculations
num_var = 3*num_nodes - 6;
num_eq = (num_nodes - 1) * num_nodes / 2;

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = field_size * rand(num_var + 6, 1);
% for i = 1:(length(x_real)/3)
%     x_real(3*i) = x_real(3*i) / 10;
% end

for loops = 1:5
    
    %Temp:
    cost = [];
    x_guess = [];
    lambda = [];
    
x_guess{1} = field_size * rand(num_var, 1);
% Temporary: guess based on real
% x_guess{1} = x_real + 10*rand(size(x_real));
% x_guess{1}([1 2 3 5 6 9]) = [];

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
        if mod(branch_per_node, 2) == 0
            neighbors = mod(i + (-branch_per_node/2:branch_per_node/2) - 1, num_nodes) + 1;
        else
            neighbors = mod(i + (-(branch_per_node-1)/2:(branch_per_node-1)/2) - 1, num_nodes) + 1;
            add_node = mod(i + num_nodes/2 - 1, num_nodes) + 1;
            neighbors = union(neighbors, add_node);
        end
        
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
lambda(1) = lambda_init;

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
    
    % Levenberg-Marquadt step
    unadjustedStep = x_guess{iter-1} - (J' * J + lambda(iter-1) * eye(size(J,2))) \ J' * s;
    adjustedStep = x_guess{iter-1} - (J' * J + lambda(iter-1) * eye(size(J,2)) / nu) \ J' * s;
    
    % Sum of squares
    [s_un, c_un] = sumOfSquares(unadjustedStep, n, m, n_ind, m_ind, null_branches, r_sq);
    [s_ad, c_ad] = sumOfSquares(adjustedStep, n, m, n_ind, m_ind, null_branches, r_sq);
    
    % Update lambda (Levenberg-Marquadt parameter)
    if c_ad < cost(iter-1)
        lambda(iter) = lambda(iter-1) / nu;
        cost(iter) = c_ad;
        s = s_ad;
        x_guess{iter} = adjustedStep;
    else
        if c_un < cost(iter-1)
            lambda(iter) = lambda(iter-1);
            cost(iter) = c_un;
            s = s_un;
            x_guess{iter} = unadjustedStep;
        else
            for k = 1:5
                lambda(iter) = lambda(iter-1) * nu^k;
                loopStep = x_guess{iter-1} - (J' * J + lambda(iter) * eye(size(J,2))) \ J' * s;
                [s_loop, c_loop] = sumOfSquares(loopStep, n, m, n_ind, m_ind, null_branches, r_sq);
                if c_loop < cost(iter-1)
                    break;
                end
            end
            cost(iter) = c_loop;
            s = s_loop;
            x_guess{iter} = loopStep;
        end
    end
    
    % Break if cost difference is low
    cost_diff = log10(abs(cost(iter) - cost(iter-1)));
    if ((cost_diff < cost_diff_thresh) && cost(iter) < cost_thresh)
        disp('Loop broken.');
        str = sprintf('Iterations: %d\nCost: %d', iter, cost(iter));
        disp(str);
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

str = sprintf('RMS error: %d', rms(e));
disp(str);

end

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

%
figure;
plot(log_cost_delta);
title('Cost Delta');
% hold on;
% plot(ones(size(cost))*thresh);
grid on;
%}

%
figure;
plot(log10(cost));
title('Cost Function');
grid on;
%}

% Lambda parameter over iterations
figure;
plot(log10(lambda));
title('Lambda per Iteration');
grid on;


%% Function
function [s, c] = sumOfSquares(guess, n, m, n_ind, m_ind, null_branches, r_sq)

    s = zeros(length(n),1);
    for i = 1:length(n)
        if (n(i) == 2 && m(i) == 1)
            s(i) = guess(1)^2 - r_sq(i);
        elseif (n(i) == 3 && m(i) == 1)
            s(i) = sum(guess(2:3).^2) - r_sq(i);
        elseif (m(i) == 1)
            s(i) = sum(guess(n_ind(i,1:3)).^2) - r_sq(i);
        elseif (n(i) == 3 && m(i) == 2)
            s(i) = (guess(2) - guess(1))^2 + guess(3)^2 - r_sq(i);
        elseif (m(i) == 2)
            s(i) = (guess(n_ind(i,1)) - guess(1))^2 + sum(guess(n_ind(i,2:3)).^2) - r_sq(i);
        elseif (m(i) == 3)
            s(i) = sum((guess(n_ind(i,1:2)) - guess(2:3)).^2) + guess(n_ind(i,3))^2 - r_sq(i);
        else
            s(i) = sum((guess(n_ind(i,:)) - guess(m_ind(i,:))).^2) - r_sq(i);
        end
    end
    
    % Nullify branches
    s(null_branches) = [];
    
    c = sum(s.^2);
end





