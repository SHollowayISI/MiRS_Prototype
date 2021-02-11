clear variables;
close all;

% Properties
num_nodes = 10;
num_iter = 5;

% Calculations
num_var = 3*num_nodes;
num_eq = (num_nodes - 1) * num_nodes / 2;

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = 1 * rand(num_var, 1);
x_guess{1} = 1 * rand(num_var, 1);
x_guess{1}([1 2 3 5 6 9]) = 0;

% Calculate exact ranges (squared)
n_ind = 3*floor((n-1)) + (1:3);
m_ind = 3*floor((m-1)) + (1:3);
r_sq = sum((x_real(n_ind) - x_real(m_ind)).^2, 2);

s = sum((x_guess{1}(n_ind) - x_guess{1}(m_ind)).^2, 2) - r_sq;
cost{1} = sum(s.^2);

for iter = 2:num_iter
    
    % Jacobian
    J = zeros(num_eq, num_var);
    for row = 1:num_eq
        J(row, n_ind(row,:)) = 2 * (x_guess{iter-1}(n_ind(row,:)) - x_guess{iter-1}(m_ind(row,:)));
        J(row, m_ind(row,:)) = 2 * (x_guess{iter-1}(m_ind(row,:)) - x_guess{iter-1}(n_ind(row,:)));
    end
    
    % Gauss-Newton step
    x_guess{iter} = x_guess{iter-1} - inv(J' * J) * J' * s;
    
    % Sum of squares
    s = sum((x_guess{iter}(n_ind) - x_guess{iter}(m_ind)).^2, 2) - r_sq;
    cost{iter} = sum(s.^2);
    
    % Print out cost
    disp(cost{iter})
    
end

% Plot points
figure;
x_plot = reshape(x_real, num_nodes, 3);
scatter3(x_plot(:,1), x_plot(:,2), x_plot(:,3), '.', 'r');
hold on;
x_plot = reshape(x_guess{end}, num_nodes, 3);
scatter3(x_plot(:,1), x_plot(:,2), x_plot(:,3));













