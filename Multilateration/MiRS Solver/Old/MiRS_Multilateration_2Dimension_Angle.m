
%% Bookkeeping

clear variables;
close all;

%% Setup

% Simulation variables
num_nodes = 10;
field_size = 100;
error_var_r = 0.0025;
error_var_a = 0.2;

% Solver variables
num_iter = 1000;
lambda_init = 100;
nu = 1.1;
r_weight = 1/(error_var_r)^0.5;
a_weight = 1/(error_var_a)^0.5;
cost_thresh = Inf;
cost_diff_thresh = -10;

% Calculations
num_var = 2*(num_nodes - 1);
num_eq = (num_nodes - 1) * num_nodes / 2;

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = field_size * (rand(num_var, 1));
x_real = [0; 0; x_real];
x_guess{1} = x_real(3:end) + 0.01*field_size * (rand(num_var, 1) - 0.5);
% x_guess{1} = field_size * (rand(num_var, 1));

% Calculate exact ranges (squared)
n_ind = 2*floor((n-1)) + (1:2) - 2;
m_ind = 2*floor((m-1)) + (1:2) - 2;
r_sq = sum((x_real(n_ind + 2) - x_real(m_ind + 2)).^2, 2);

% Calculate exact angles
ang = atan2d(x_real(n_ind(:,2) + 2) - x_real(m_ind(:,2) + 2), x_real(n_ind(:,1) + 2) - x_real(m_ind(:,1) + 2));

% Introduce error into range and angle measurement
r_sq = r_sq + error_var_r*(rand(size(r_sq)) - 0.5);
ang = ang + error_var_a*(rand(size(ang)) - 0.5);

%% Initial step

% Evaluate cost function
[s, c, r] = sumOfSquares(x_guess{1}, n, m, n_ind, m_ind, r_sq, ang);

cost(1) = c;
lambda(1) = lambda_init;

%% Main loop

for iter = 2:num_iter
    
    % Jacobian
    J = zeros(2*num_eq, num_var);
    for row = 1:num_eq
        if m(row) == 1
            J(row, n_ind(row, :)) = 2 * (x_guess{iter-1}(n_ind(row,:)));
            J(row + num_eq, n_ind(row,1)) = -x_guess{iter-1}(n_ind(row,2));
            J(row + num_eq, n_ind(row,2)) = x_guess{iter-1}(n_ind(row,1));
        else
            J(row, n_ind(row,:)) = 2 * (x_guess{iter-1}(n_ind(row,:)) - x_guess{iter-1}(m_ind(row,:)));
            J(row, m_ind(row,:)) = 2 * (x_guess{iter-1}(m_ind(row,:)) - x_guess{iter-1}(n_ind(row,:)));
            J(row + num_eq, n_ind(row,1)) = x_guess{iter-1}(m_ind(row,2)) - x_guess{iter-1}(n_ind(row,2));
            J(row + num_eq, m_ind(row,1)) = x_guess{iter-1}(n_ind(row,2)) - x_guess{iter-1}(m_ind(row,2));
            J(row + num_eq, n_ind(row,2)) = x_guess{iter-1}(n_ind(row,1)) - x_guess{iter-1}(m_ind(row,1));
            J(row + num_eq, m_ind(row,2)) = x_guess{iter-1}(m_ind(row,1)) - x_guess{iter-1}(n_ind(row,1));
        end
    end
    J((num_eq+1):end,:) = rad2deg(J((num_eq+1):end,:) ./ r);
    
    % Create weights
    W = diag([r_weight*ones(num_eq, 1); a_weight*ones(num_eq, 1)]);
    
%     x_guess{iter} = x_guess{iter-1} - (J' * W * J + lambda_init * diag(J' * W * J) .* eye(size(J,2))) \ J' * W * s;
%     [s, cost(iter), r] = sumOfSquares(x_guess{iter}, n, m, n_ind, m_ind, r_sq, ang);
    
    % Levenberg-Marquadt step
    unadjustedStep = x_guess{iter-1} + (J' * W * J + lambda(iter-1) * diag(J' * W * J) .* eye(size(J,2))) \ J' * W * s;
    adjustedStep = x_guess{iter-1} + (J' * W * J + lambda(iter-1) * diag(J' * W * J) .* eye(size(J,2)) / nu) \ J' * W * s;
        
    % Sum of squares
    [s_un, c_un, r_un] = sumOfSquares(unadjustedStep, n, m, n_ind, m_ind, r_sq, ang);
    [s_ad, c_ad, r_ad] = sumOfSquares(adjustedStep, n, m, n_ind, m_ind, r_sq, ang);
    
    % Update lambda (Levenberg-Marquadt parameter)
    if c_ad < cost(iter-1)
        lambda(iter) = lambda(iter-1) / nu;
        cost(iter) = c_ad;
        s = s_ad;
        r = r_ad;
        x_guess{iter} = adjustedStep;
    else
        if c_un < cost(iter-1)
            lambda(iter) = lambda(iter-1);
            cost(iter) = c_un;
            s = s_un;
            r = r_un;
            x_guess{iter} = unadjustedStep;
        else
            for k = 1:5
                lambda(iter) = lambda(iter-1) * nu^k;
                loopStep = x_guess{iter-1} + ((J' * W * J + lambda(iter) * diag(J' * W * J) .* eye(size(J,2)))) \ J' * W * s;
                [s_loop, c_loop, r_loop] = sumOfSquares(loopStep, n, m, n_ind, m_ind, r_sq, ang);
                if c_loop < cost(iter-1)
                    break;
                end
            end
            cost(iter) = c_loop;
            s = s_loop;
            r = r_loop;
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
x_result = reshape([0; 0; x_guess{end}], 2, num_nodes)';
x_exact = reshape(x_real, 2, num_nodes)';

% Perform rotation
% x1 = x_exact(2,1);
% y1 = x_exact(2,2);
% x2 = x_result(2,1);
% y2 = x_result(2,2);
% 
% theta = atan2d(y2, x2) - atan2d(y1, x1);
theta = mean(atan2d(x_result(2:end,2), x_result(2:end,1)) - atan2d(x_exact(2:end,2), x_exact(2:end,1)));

R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
x_rotated = transpose(R * x_exact');

% Calculate error
e = sum((x_rotated - x_result).^2, 2);
str = sprintf('\nAngle Error: %0.3f [deg]\nRMS Distance Error: %d [m], \nMax Distance Error: %d [m]', theta, rms(e), max(e));
disp(str);

% Plotting
close all;

% Plot points
figure;
scatter(x_result(:,1), x_result(:,2), '.', 'r');
hold on;
scatter(x_exact(:,1), x_exact(:,2));
hold on;
scatter(x_rotated(:,1), x_rotated(:,2), 'b');
grid on;
xlim([0 1] * field_size)
ylim([0 1] * field_size)


%% Cost function
function [s, c, r] = sumOfSquares(x_in, n, m, n_ind, m_ind, r_sq, ang)

% Evaluate cost function
r_sq_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        r_sq_est(i) = sum(x_in(n_ind(i,:)).^2);
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

s = [ang_est - ang; r_sq_est - r_sq];

r = r_sq_est;

c = sum(s.^2);

end












