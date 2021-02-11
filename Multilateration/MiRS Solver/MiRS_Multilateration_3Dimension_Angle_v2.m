
%% Bookkeeping

clear variables;
close all;

%% User Inputs

% Simulation variables
num_nodes = 10;
field_size_xy = 1000;
field_size_z = 10;
offset_z = 0;
error_var_r = 0.0025;
error_var_a = 0.2;

% Solver variables
num_iter = 1000;
r_factor = 1e16;
a_factor = 1e-4;
m_factor = 1e16;
c_factor = 1e12;
z_factor = 1e12;
r_weight = r_factor/(error_var_r)^4;
a_weight = a_factor/(error_var_a)^2;
z_weight = z_factor/(field_size_z)^2;
c_weight = c_factor * 12/((num_nodes - 1)*(field_size_z^2));
cost_thresh = 0;
cost_diff_thresh = -Inf;
lambda_init = 0.01;
L_up = 11;
L_down = 9;
epsilon = 0.9;

%% Setup

% Calculations
num_var = 3*(num_nodes - 1);
num_eq = (num_nodes - 1) * num_nodes / 2;

% Initialize loop variables
cost = zeros(num_iter, 1);
lambda = zeros(num_iter, 1);
x_guess = repmat({zeros(num_var, 1)}, num_iter, 1);

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = field_size_xy * (rand(num_var, 1) - 0.5);
x_real(mod(1:length(x_real), 3) == 0) = field_size_z * (rand((num_nodes-1),1) - 0.5);
x_real = [0; 0; offset_z; x_real];

% Calculate exact ranges (squared)
n_ind = 3*floor((n-1)) + (1:3) - 3;
m_ind = 3*floor((m-1)) + (1:3) - 3;

% Calculate exact angles and ranges
r_sq = sum((x_real(n_ind + 3) - x_real(m_ind + 3)).^2, 2);
ang = atan2d(x_real(n_ind(:,2) + 3) - x_real(m_ind(:,2) + 3), x_real(n_ind(:,1) + 3) - x_real(m_ind(:,1) + 3));

% First guess using blurry real points from master node
guess_ang = ang(m == 1) + 2 * error_var_a * randn(size(ang(m == 1)));
guess_r = sqrt(r_sq(m == 1)) + 2 * error_var_r * randn(size(r_sq(m == 1)));
guess_x = guess_r .* cosd(guess_ang);
guess_y = guess_r .* sind(guess_ang);
guess_z = field_size_z * (rand(size(guess_ang))-0.5);
x_guess{1} = reshape(([guess_x, guess_y, guess_z])', [], 1);

% Introduce error into range and angle measurement
r_sq = r_sq + error_var_r*(rand(size(r_sq)) - 0.5);
ang = ang + error_var_a*(rand(size(ang)) - 0.5);
z_s = guess_z;
c_s = 0;

% Vectorize data
% y = [r_sq; ang];
% y = [r_sq; ang; z_s];
% y = [r_sq; ang; c_s];
y = [r_sq; ang; z_s; c_s];

%% Initial step

% Matrix setup
m_weight = ones(num_eq, 1);
m_weight(m == 1) = m_factor;
% W = diag([r_weight*ones(num_eq, 1); a_weight*m_weight]);
% W = diag([r_weight*ones(num_eq, 1); a_weight*m_weight; z_weight*ones(num_var/3, 1)]);
% W = diag([r_weight*ones(num_eq, 1); a_weight*m_weight; c_weight]);
W = diag([r_weight*ones(num_eq, 1); a_weight*m_weight; z_weight*ones(num_var/3, 1); c_weight]);
J = jacobian(x_guess{1}, n, m, n_ind, m_ind, num_var, num_eq);

% Evaluate cost function
y_hat = angleRangeFunctions(x_guess{1}, n, m, n_ind, m_ind, offset_z);
cost(1) = costFunction(y_hat, y, W, J, zeros(num_var, 1));

% Initialize lambda parameter
lambda(1) = lambda_init;

%% Main loop

for iter = 2:num_iter
    
    % Get Jacobian matrix
    J = jacobian(x_guess{iter-1}, n, m, n_ind, m_ind, num_var, num_eq);
    
    % Levenberg-Marquadt Step
    h = LMStep(y_hat, y, W, J, lambda(iter-1));
    
    % Lambda iteration
%     eta(iter) = LMMetric(y_hat, y, W, J, h, lambda(iter-1));
    if LMMetric(y_hat, y, W, J, h, lambda(iter-1)) > epsilon
        cost(iter) = costFunction(y_hat, y, W, J, h);
        x_guess{iter} = x_guess{iter-1} + h;
        lambda(iter) = max(lambda(iter-1)/L_down, 1e-7);
    else
        cost(iter) = cost(iter-1);
        x_guess{iter} = x_guess{iter-1};
        lambda(iter) = min(lambda(iter-1)*L_up, 1e7);
    end
    
    % Update guess and cost
    cost(iter) = costFunction(y_hat, y, W, J, h);
    x_guess{iter} = x_guess{iter-1} + h;
    y_hat = angleRangeFunctions(x_guess{iter}, n, m, n_ind, m_ind, offset_z);
    
    % Temp: retain lambda
    lambda(iter) = lambda(iter-1);
    
    % Break if cost difference is low
    cost_diff = log10(abs(cost(iter) - cost(iter-1)));
    if ((cost_diff < cost_diff_thresh) && (cost(iter) < cost_thresh))
        disp('Loop broken.');
        fprintf('Iterations: %d\nCost: %d\n', iter, cost(iter));
        break;
    end
    
end

%% Error estimation & Plotting

% Degree of freedom reduction and error estimation
x_result = reshape([0; 0; offset_z; x_guess{end}], 3, num_nodes)';
x_exact = reshape(x_real, 3, num_nodes)';

% Updated rotation method
com = mean(x_result(2:end,:), 1);
x = x_result(2:end,1) - com(1);
y = x_result(2:end,2) - com(2);
z = x_result(2:end,3) - com(3);
D = (x'*x)*(y'*y) - (x'*y)*(x'*y);
a = ((y'*z)*(x'*y) - (x'*z)*(y'*y))/D;
b = ((x'*y)*(x'*z) - (x'*x)*(y'*z))/D;
n = [a, b, 1];
theta = -acosd(1/sqrt(sum(n.^2)));
k = [-b, a, 0];
k = k ./ norm(k);
K = [0, 0, k(2); 0, 0, -k(1); -k(2), k(1), 0];
R = eye(3) + sind(theta)*K + (1-cosd(theta))*K*K;
x_rotated = zeros(size(x_result));
x_rotated(1,:) = x_result(1,:);
x_rotated(2:end,:) = transpose(R * x_result(2:end,:)');

% Perform rotation
phi = -mean(atan2d(x_rotated(2:end,2), x_rotated(2:end,1)) - atan2d(x_exact(2:end,2), x_exact(2:end,1)));
R_z = [cosd(phi), -sind(phi), 0; sind(phi), cosd(phi), 0; 0, 0, 1];
x_rotated = transpose(R_z * x_rotated');

% Calculate error
error_list = x_rotated(2:end,:) - x_exact(2:end,:);
err_off = sqrt(sum(mean(error_list,1).^2));
err_rms = rms(sqrt(sum((error_list - mean(error_list,1)).^2, 2)));
err_xy_off = sqrt(sum(mean(error_list(:,1:2),1).^2));
err_xy_rms = rms(sqrt(sum((error_list(:,1:2) - mean(error_list(:,1:2),1)).^2, 2)));
err_z_off = mean(error_list(:,3),1);
err_z_rms = rms(error_list(:,3) - mean(error_list(:,3),1));

fprintf('\nSimulation Complete.\n');
fprintf('Angle Error: (%0.3f x %0.3f) [deg]\n', theta, phi);
fprintf('\nTotal Distance: \nError Offset: %d [m], \nError RMS: %d [m]\n', err_off, err_rms);
fprintf('\nXY-Plane: \nError Offset: %d [m], \nError RMS: %d [m]\n', err_xy_off, err_xy_rms);
fprintf('\nAltitude: \nError Offset: %d [m], \nError RMS: %d [m]\n', err_z_off, err_z_rms);


% Plotting
close all;

% Plot points
figure;
scatter3(x_result(:,1), x_result(:,2), x_result(:,3), '.', 'r');
hold on;
scatter3(x_exact(:,1), x_exact(:,2), x_exact(:,3));
hold on;
scatter3(x_rotated(:,1), x_rotated(:,2), x_rotated(:,3), '.', 'b');
% hold on;
% scatter3(x_rotated_2(:,1), x_rotated_2(:,2), x_rotated_2(:,3), 'g');
grid on;
xlim([-1 1] * field_size_xy)
xlabel('x')
ylim([-1 1] * field_size_xy)
ylabel('y')
% zlim([-field_size_z field_size_z + offset_z])
zlim([-1 1] * field_size_xy)

%% Cost function
function y_hat = angleRangeFunctions(x_in, n, m, n_ind, m_ind, offset_z)

% Evaluate cost function
r_sq_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        r_sq_est(i) = sum(x_in(n_ind(i,1:2)).^2) + (x_in(n_ind(i,3)) - offset_z).^2;
    else
        r_sq_est(i) = sum((x_in(n_ind(i,:)) - x_in(m_ind(i,:))).^2);
    end
end

% Cost function for angle
ang_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        ang_est(i) = atan2d(x_in(n_ind(i,2)), x_in(n_ind(i,1)));
    else
        ang_est(i) = atan2d(x_in(n_ind(i,2)) - x_in(m_ind(i,2)), x_in(n_ind(i,1)) - x_in(m_ind(i,1)));
    end
end

% % Cost function for altitude
z_est = zeros(length(x_in)/3, 1);
for i = 1:(length(x_in)/3)
    z_est(i) = x_in(3*i);
end

% Cost function for z-center
c_est = (sum(x_in(mod(1:length(x_in),3) == 0)) / (length(x_in)/3)).^2;

% y_hat = [r_sq_est; ang_est];
% y_hat = [r_sq_est; ang_est; z_est];
% y_hat = [r_sq_est; ang_est; c_est];
y_hat = [r_sq_est; ang_est; z_est; c_est];

end

%% Jacobian
function J = jacobian(x_in, n, m, n_ind, m_ind, num_var, num_eq)

% Calculate squared range between points
r_sq_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        r_sq_est(i) = sum(x_in(n_ind(i,:)).^2);
    else
        r_sq_est(i) = sum((x_in(n_ind(i,:)) - x_in(m_ind(i,:))).^2);
    end
end

% Calculate jacobian
% J = zeros(2*num_eq, num_var);
% J = zeros(2*num_eq + (num_var/3), num_var);
% J = zeros(2*num_eq + 1, num_var);
J = zeros(2*num_eq + (num_var/3) + 1, num_var);
for row = 1:num_eq
    if m(row) == 1
        J(row, n_ind(row, :)) = 2 * (x_in(n_ind(row,:)));
        J(row + num_eq, n_ind(row,1)) = -x_in(n_ind(row,2));
        J(row + num_eq, n_ind(row,2)) = x_in(n_ind(row,1));
    else
        J(row, n_ind(row,:)) = 2 * (x_in(n_ind(row,:)) - x_in(m_ind(row,:)));
        J(row, m_ind(row,:)) = 2 * (x_in(m_ind(row,:)) - x_in(n_ind(row,:)));
        J(row + num_eq, n_ind(row,1)) = x_in(m_ind(row,2)) - x_in(n_ind(row,2));
        J(row + num_eq, m_ind(row,1)) = x_in(n_ind(row,2)) - x_in(m_ind(row,2));
        J(row + num_eq, n_ind(row,2)) = x_in(n_ind(row,1)) - x_in(m_ind(row,1));
        J(row + num_eq, m_ind(row,2)) = x_in(m_ind(row,1)) - x_in(n_ind(row,1));
    end
end
J((num_eq+1):(2*num_eq),:) = rad2deg(J((num_eq+1):(2*num_eq),:) ./ r_sq_est);

% Z constraint
for row = 1:(num_var/3)
    J(row + 2*num_eq, 3*row) = 1;
end

% Z-center constraint
for row = 1:(num_var/3)
    J(end, 3*row) = 2 * sum(x_in(mod(1:num_var, 3) == 0)) / ((num_var/3)^2);
end

end

%% Cost function
function chi_sq = costFunction(y_hat, y, W, J, h)

% Calculate chi-squared cost function value
chi_sq = (y' * W * y) + (y_hat' * W * y_hat) - 2 * (y' * W * y_hat) + ...
    (h' * J' * W * J * h) - 2 * ((y - y_hat)' * W * J * h);

end

%% Levenberg-Marquadt Step
function h = LMStep(y_hat, y, W, J, lambda)

% Calculate step
h = ((J' * W * J) + lambda * diag(J' * W * J) .* eye(size(J, 2))) \ (J' * W * (y - y_hat));


end

%% Lambda Update Metric
function rho = LMMetric(y_hat, y, W, J, h, lambda)

% Calculate metric for L-M parameter update
rho = (costFunction(y_hat, y, W, J, zeros(size(h))) - costFunction(y_hat, y, W, J, h)) ./ ...
    (h' * (lambda * diag(J' * W * J) .* eye(size(J, 2)) * h + (J' * W * (y - y_hat))));

end









