%% SF2527 Black--Scholes project
clear;
clc;
close all;

%% Parameters from the paper
K = 100;
r = 0.03;
sigma = 0.3;
T = 0.1;
S0 = 100;
h = 0.05;
Smax = 400;
Nt_values = [160, 320, 640, 1280];
plot_limit = 160;

%% Article convergence table
exact_price = black_scholes_call(S0, T, K, r, sigma);
paper_price = [3.1676, 3.7501, 3.8871, 3.9189];
ade_price = zeros(size(Nt_values));
exact_error = zeros(size(Nt_values));

for q = 1:length(Nt_values)
    Nt = Nt_values(q);
    [x, u, tau_q, history_q] = ...
        ade_method(K, r, sigma, T, h, Smax, Nt, plot_limit);

    [~, target_index] = min(abs(x - S0));
    ade_price(q) = u(target_index);
    exact_error(q) = abs(ade_price(q) - exact_price);

    if q == length(Nt_values)
        x_ade = x;
        u_ade = u;
        tau = tau_q;
        u_ade_history = history_q;
    end
end

order = [NaN, log2(exact_error(1:end-1)./exact_error(2:end))];
results = table(Nt_values.', ade_price.', paper_price.', exact_error.', order.', ...
    'VariableNames', {'Nt', 'ADEPrice', 'PaperPrice', 'ExactError', 'Order'});
disp(results)

%% Space convergence
h_values = [0.8, 0.4, 0.2, 0.1];
Nt_space = [80, 320, 1280, 5120]; % dt is proportional to h^2
space_error = zeros(size(h_values));

for q = 1:length(h_values)
    [x_space, u_space] = ade_method( ...
        K, r, sigma, T, h_values(q), 800, Nt_space(q), plot_limit);
    [~, target_index] = min(abs(x_space - S0));
    space_error(q) = abs(u_space(target_index) - exact_price);
end

space_order = [NaN, log2(space_error(1:end-1)./space_error(2:end))];
space_results = table(h_values.', Nt_space.', space_error.', space_order.', ...
    'VariableNames', {'h', 'Nt', 'Error', 'Order'});
disp(space_results)

figure;
loglog(h_values, space_error, 'o-', 'LineWidth', 1.3)
xlabel('Space step h')
ylabel('Error at S = 100')
title('Space convergence')
grid on

%% Forward Euler stability comparison
[x_euler, u_euler, tau_euler, ~] = ...
    forward_euler(K, r, sigma, T, h, Smax, Nt_values(end), plot_limit);
fprintf('Forward Euler stopped at tau = %.6f before T = %.1f because it became unstable.\n', ...
    tau_euler(end), T)

%% Compare with the exact solution
show_ade = x_ade <= plot_limit;
exact_profile = black_scholes_call(x_ade(show_ade), T, K, r, sigma);

figure;
plot(x_ade(show_ade), exact_profile, 'k-', 'LineWidth', 1.5)
hold on
plot(x_ade(show_ade), u_ade(show_ade), '--', 'LineWidth', 1.3)
hold off
xlabel('Asset price')
ylabel('Call option value')
legend('Exact', 'ADE method', 'Location', 'northwest')
grid on

%% ADE solution over asset price and time
x_ade_plot = x_ade(1:size(u_ade_history, 2));

figure;
mesh(x_ade_plot(1:20:end), tau(1:16:end), u_ade_history(1:16:end, 1:20:end))
title('ADE method')
xlabel('Asset price')
ylabel('Time to maturity')
zlabel('Call option value')
xlim([0, plot_limit])
grid on

%% Forward Euler instability
show_euler = x_euler <= plot_limit;
exact_euler = black_scholes_call(x_euler(show_euler), tau_euler(end), K, r, sigma);

figure;
plot(x_euler(show_euler), exact_euler, 'k-', 'LineWidth', 1.5)
hold on
plot(x_euler(show_euler), u_euler(show_euler), '--', 'LineWidth', 1.3)
hold off
title(sprintf('Forward Euler at tau = %.4f', tau_euler(end)))
xlabel('Asset price')
ylabel('Call option value')
legend('Exact', 'Forward Euler', 'Location', 'northwest')
ylim([-20, 80])
grid on

function price = black_scholes_call(S, tau, K, r, sigma)
d1 = (log(S/K) + (r + sigma^2/2)*tau)/(sigma*sqrt(tau));
d2 = d1 - sigma*sqrt(tau);
Phi = @(z) 0.5*erfc(-z/sqrt(2));
price = S.*Phi(d1) - K*exp(-r*tau).*Phi(d2);
end
