
%% Bookkeeping

clear variables;
close all;

%% Setup

% Properties
num_nodes = 20;
num_iter = 100;
field_size = 100;
error_var = 1;

% Calculations
num_var = 2*num_nodes - 3;
num_eq = (num_nodes - 1) * num_nodes / 2;

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = field_size * rand(num_var + 3, 1);
x_guess{1} = field_size * rand(num_var + 3, 1);

% Remove fixed variables
x_guess{1}(1:3) = [];

% Calculate exact ranges (squared)
n_ind = 2*floor((n-1)) + (1:2) - 3;
m_ind = 2*floor((m-1)) + (1:2) - 3;
r_sq = sum((x_real(n_ind + 3) - x_real(m_ind + 3)).^2, 2);

% Introduce error into range measurement
r_sq = r_sq + error_var*rand(size(r_sq));

%% Initial step

% Evaluate cost function
s = zeros(length(n),1);
for i = 1:length(n)
    if (n(i) == 2)
        s(i) = x_guess{1}(n_ind(i,2)).^2 - r_sq(i);
    elseif (m(i) == 1)
        s(i) = sum(x_guess{1}(n_ind(i,:)).^2) - r_sq(i);
    elseif (m(i) == 2)
        s(i) = x_guess{1}(n_ind(i,1))^2 +  (x_guess{1}(n_ind(i,2)) - x_guess{1}(m_ind(i,2)))^2 - r_sq(i); 
    else
        s(i) = sum((x_guess{1}(n_ind(i,:)) - x_guess{1}(m_ind(i,:))).^2) - r_sq(i);
    end
end

cost{1} = sum(s.^2);

%% Main loop

for iter = 2:num_iter
    
    % Jacobian
    J = zeros(num_eq, num_var);
    for row = 1:num_eq
        if n(row) == 2
            J(row, n_ind(row,2)) = 2 * (x_guess{iter-1}(n_ind(row,2)));
        elseif m(row) == 1
            J(row, n_ind(row, :)) = 2* (x_guess{iter-1}(n_ind(row,:)));
        elseif m(row) == 2
            J(row, n_ind(row, 1)) = 2* (x_guess{iter-1}(n_ind(row,1)));
            J(row, n_ind(row, 2)) = 2* (x_guess{iter-1}(n_ind(row,2)) - x_guess{iter-1}(m_ind(row,2)));
            J(row, m_ind(row, 2)) = 2* (x_guess{iter-1}(m_ind(row,2)) - x_guess{iter-1}(n_ind(row,2)));
        else
            J(row, n_ind(row,:)) = 2 * (x_guess{iter-1}(n_ind(row,:)) - x_guess{iter-1}(m_ind(row,:)));
            J(row, m_ind(row,:)) = 2 * (x_guess{iter-1}(m_ind(row,:)) - x_guess{iter-1}(n_ind(row,:)));
        end
    end
    
    % Gauss-Newton step
    x_guess{iter} = x_guess{iter-1} - inv(J' * J) * J' * s;
    
    % Sum of squares
    s = zeros(length(n),1);
    for i = 1:length(n)
        if (n(i) == 2)
            s(i) = x_guess{iter}(n_ind(i,2)).^2 - r_sq(i);
        elseif (m(i) == 1)
            s(i) = sum(x_guess{iter}(n_ind(i,:)).^2) - r_sq(i);
        elseif (m(i) == 2)
            s(i) = x_guess{iter}(n_ind(i,1))^2 +  (x_guess{iter}(n_ind(i,2)) - x_guess{iter}(m_ind(i,2)))^2 - r_sq(i);
        else
            s(i) = sum((x_guess{iter}(n_ind(i,:)) - x_guess{iter}(m_ind(i,:))).^2) - r_sq(i);
        end
    end
    cost{iter} = sum(s.^2);
    
end

%% Error estimation & Plotting

% Degree of freedom reduction and error estimation
x_result = reshape([0; 0; 0; x_guess{end}], 2, num_nodes)';
x_exact = reshape(x_real, 2, num_nodes)';
x_exact = x_exact - x_exact(1,:);

% Perform rotation
x1 = x_exact(2,1);
y1 = x_exact(2,2);
y2 = x_result(2,2);

theta = atan2d(x1*y2, y1*y2);
R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
x_rotated = transpose(R * x_exact');

% Perform reflection
flip = sign(x_rotated(3,1))/sign(x_result(3,1));
x_rotated = x_rotated .* [flip, 1];

% Calculate error
e = sum((x_rotated - x_result).^2, 2);

rms(e)

% Plotting
close all;

% Plot points
figure;
scatter(x_result(:,1), x_result(:,2), '.', 'r');
hold on;
scatter(x_rotated(:,1), x_rotated(:,2));
grid on;
xlim([-1 1] * field_size)
ylim([-1 1] * field_size)













