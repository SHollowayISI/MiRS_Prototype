
x = x_result(:,1);
y = x_result(:,2);
z = x_result(:,3);

V = [x.*z; y.*z; z];
M = [x.*x, x.*y, x; x.*y, y.*y, y; x, y, ones(size(x))];

A = M \ V;

n = [A(1); A(2); 1];

k = cross(n, [0; 0; 1]);
k = k / norm(k);

theta = acosd(1 / norm(n));

k_mat = [0, 0, k(2); 0, 0, -k(1); -k(2), k(1), 0];

R = eye(3) + sind(theta) * k_mat + (1-cosd(theta)) * k_mat * k_mat;

x_rotated = transpose(R * x_result');

% Plotting
close all;

% Plot points
figure;
scatter3(x_result(:,1), x_result(:,2), x_result(:,3), '.', 'r');
hold on;
scatter3(x_exact(:,1), x_exact(:,2), x_exact(:,3));
hold on;
scatter3(x_rotated(:,1), x_rotated(:,2), x_rotated(:,3), 'b');
grid on;
xlim([-1 1] * field_size_xy)
xlabel('x')
ylim([-1 1] * field_size_xy * 0.8)
ylabel('y')
zlim([-field_size_z-offset_z field_size_z + offset_z])