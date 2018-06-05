close all
clear all
tic
num_runners = 200;
x = unifrnd(-1, 0, num_runners, 1);

sd = 0.5
mu = 5.
v_max = sqrt(sd) * randn(size(x)) + mu;
v_min = 1.


rhs = @(t, x) v_min + v_max .* max((1-density(x)./10.), 0.);


[ts, y] = ode45(rhs, [0 100], x);
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

