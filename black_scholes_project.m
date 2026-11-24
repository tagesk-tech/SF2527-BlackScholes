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

%% Functions
function [x, u, tau, u_history] = ade_method(K, r, sigma, T, h, Smax, Nt, plot_limit)
x = 0:h:Smax;
u = max(x - K, 0);
phi = zeros(size(u));
psi = zeros(size(u));

dt = T/Nt;
m = length(x);
plot_m = min(find(x <= plot_limit, 1, 'last'), m - Nt);
tau = (0:Nt)*dt;
u_history = zeros(Nt + 1, plot_m);
u_history(1, :) = u(1:plot_m);

for n = 1:Nt
    mnew = m - 1;
    phi(1) = 0;

    % Left-to-right sweep, equation (8)
    for i = 2:mnew
        diffusion = sigma^2*x(i)^2/(2*h^2);
        advection = r*x(i)/(2*h);

        numerator = u(i)/dt ...
            + diffusion*(u(i+1) - u(i) + phi(i-1)) ...
            + advection*(u(i+1) - u(i) - phi(i-1)) ...
            - r*u(i)/2;
        denominator = 1/dt + diffusion - advection + r/2;
        phi(i) = numerator/denominator;
    end

    % Right-to-left sweep, equation (9)
    psi(mnew) = phi(mnew);
    for i = mnew-1:-1:2
        diffusion = sigma^2*x(i)^2/(2*h^2);
        advection = r*x(i)/(2*h);

        numerator = u(i)/dt ...
            + diffusion*(psi(i+1) - u(i) + u(i-1)) ...
            + advection*(psi(i+1) + u(i) - u(i-1)) ...
            - r*u(i)/2;
        denominator = 1/dt + diffusion + advection + r/2;
        psi(i) = numerator/denominator;
    end

    psi(1) = 0;
    u(1:mnew) = (phi(1:mnew) + psi(1:mnew))/2; % equation (10)
    u_history(n + 1, :) = u(1:plot_m);
    m = mnew;
end

x = x(1:m);
u = u(1:m);
end

function [x, u, tau, u_history] = forward_euler(K, r, sigma, T, h, Smax, Nt, plot_limit)
x = 0:h:Smax;
u = max(x - K, 0);

dt = T/Nt;
m = length(x);
plot_m = min(find(x <= plot_limit, 1, 'last'), m - Nt);
tau_all = (0:Nt)*dt;
u_history = zeros(Nt + 1, plot_m);
u_history(1, :) = u(1:plot_m);

for n = 1:Nt
    mnew = m - 1;
    i = 2:mnew;
    diffusion = sigma^2*x(i).^2/(2*h^2);
    advection = r*x(i)/(2*h);

    unew = zeros(1, mnew);
    unew(i) = u(i) + dt*( ...
        diffusion.*(u(i+1) - 2*u(i) + u(i-1)) ...
        + advection.*(u(i+1) - u(i-1)) - r*u(i));

    u(1:mnew) = unew;
    u_history(n + 1, :) = u(1:plot_m);
    m = mnew;

    % Stop once the instability is visible instead of overflowing.
    if any(~isfinite(u(1:plot_m))) || max(abs(u(1:plot_m))) > 1e3
        break
    end
end

x = x(1:m);
u = u(1:m);
tau = tau_all(1:n+1);
u_history = u_history(1:n+1, :);
end

function price = black_scholes_call(S, tau, K, r, sigma)
d1 = (log(S/K) + (r + sigma^2/2)*tau)/(sigma*sqrt(tau));
d2 = d1 - sigma*sqrt(tau);
Phi = @(z) 0.5*erfc(-z/sqrt(2));
price = S.*Phi(d1) - K*exp(-r*tau).*Phi(d2);
end
