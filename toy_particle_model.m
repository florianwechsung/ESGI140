close all
clear all
tic
num_runners = 200; % people
width_0 = 7.5; % meter
start_box_density = 4; % (people/square_meter)
start_box_length = num_runners/(width_0 * start_box_density);
x = unifrnd(-start_box_length, 0, num_runners, 1);

sd = 0.5;
mu = 5.;
v_max = sqrt(sd) * randn(size(x)) + mu;
v_min = 1;

rho_2 = 2.;
rho_1 = .5;
%rhs = @(t, x) v_min + v_max .* max((1-density(x)./10.), 0.);
rhs = @(t, x) min(max(v_max + ((v_min-v_max)/(rho_2-rho_1)) .* (density(x) - rho_1), v_min), v_max)


[ts, y] = ode45(rhs, [0 30], x);
toc
hold on
for i=1:size(y,2)
    plot(ts, y(:, i))
end

threshold = 100.;

time_to_threshold = zeros(num_runners, 1);
for i=1:size(y, 2)
    time_to_threshold(i) = min(find(y(:, i) >= threshold));
end
figure()
hist(time_to_threshold)

