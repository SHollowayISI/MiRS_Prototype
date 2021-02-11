clear variables;
close all;

% Properties
num_nodes = 10;
num_iter = 20;

% Calculations
num_var = num_nodes - 1;
num_eq = (num_nodes - 1) * num_nodes / 2;

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = 1 * rand(num_var + 1, 1);
x_guess{1} = 1 * rand(num_var + 1, 1);

% Remove fixed variables
x_guess{1}(1) = [];

% Calculate exact ranges (squared)
r_sq = sum((x_real(n) - x_real(m)).^2, 2);

% Evaluate cost function
s = zeros(length(n),1);
for i = 1:length(n)
    if m(i) == 1
        s(i) = x_guess{1}(n(i)-1).^2 - r_sq(i);
    else
        s(i) = sum((x_guess{1}(n(i)-1) - x_guess{1}(m(i)-1)).^2, 2) - r_sq(i);
    end
end

cost{1} = sum(s.^2);

for iter = 2:num_iter
    
    % Jacobian
    J = zeros(num_eq, num_var);
    for row = 1:num_eq
        if m(row) == 1
            J(row, n(row,:)-1) = 2 * (x_guess{iter-1}(n(row,:)-1));
        else
            J(row, n(row,:)-1) = 2 * (x_guess{iter-1}(n(row,:)-1) - x_guess{iter-1}(m(row,:)-1));
            J(row, m(row,:)-1) = 2 * (x_guess{iter-1}(m(row,:)-1) - x_guess{iter-1}(n(row,:)-1));
        end
    end
    
    % Gauss-Newton step
    x_guess{iter} = x_guess{iter-1} - inv(J' * J) * J' * s;
    
    % Sum of squares
    s = zeros(length(n),1);
    for i = 1:length(n)
        if m(i) == 1
            s(i) = x_guess{iter}(n(i)-1).^2 - r_sq(i);
        else
            s(i) = sum((x_guess{iter}(n(i)-1) - x_guess{iter}(m(i)-1)).^2, 2) - r_sq(i);
        end
    end
    cost{iter} = sum(s.^2);
    
    % Print out cost
    disp(cost{iter})
    
end

% % Plot points
% figure;
% x_plot = reshape(x_real, num_nodes, 3);
% scatter3(x_plot(:,1), x_plot(:,2), x_plot(:,3), '.', 'r');
% hold on;
% x_plot = reshape(x_guess{end}, num_nodes, 3);
% scatter3(x_plot(:,1), x_plot(:,2), x_plot(:,3));













