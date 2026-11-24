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
Nt = 1280;
plot_limit = 160;

%% Article method
[x_ade, u_ade, tau, u_ade_history] = ...
    ade_method(K, r, sigma, T, h, Smax, Nt, plot_limit);
[x_cn, u_cn, ~, u_cn_history] = ...
    crank_nicolson(K, r, sigma, T, h, Smax, Nt, plot_limit);

[~, ade_index] = min(abs(x_ade - S0));
[~, cn_index] = min(abs(x_cn - S0));
ade_price = u_ade(ade_index);
cn_price = u_cn(cn_index);
exact_price = black_scholes_call(S0, T, K, r, sigma);
fprintf('ADE price:              %.6f\n', ade_price)
fprintf('Crank-Nicolson price:   %.6f\n', cn_price)
fprintf('Exact price:            %.6f\n', exact_price)

%% Compare with the exact solution
show_ade = x_ade <= plot_limit;
show_cn = x_cn <= plot_limit;
exact_profile = black_scholes_call(x_ade(show_ade), T, K, r, sigma);

figure;
plot(x_ade(show_ade), exact_profile, 'k-', 'LineWidth', 1.5)
hold on
plot(x_ade(show_ade), u_ade(show_ade), '--', 'LineWidth', 1.3)
plot(x_cn(show_cn), u_cn(show_cn), ':', 'LineWidth', 1.5)
hold off
xlabel('Asset price')
ylabel('Call option value')
legend('Exact', 'ADE method', 'Crank-Nicolson', 'Location', 'northwest')
grid on

%% Numerical solution over asset price and time
x_ade_plot = x_ade(1:size(u_ade_history, 2));
x_cn_plot = x_cn(1:size(u_cn_history, 2));

figure;
tiledlayout(1, 2)

nexttile
mesh(x_ade_plot(1:20:end), tau(1:16:end), u_ade_history(1:16:end, 1:20:end))
title('ADE method')
xlabel('Asset price')
ylabel('Time to maturity')
zlabel('Call option value')
xlim([0, plot_limit])
grid on

nexttile
mesh(x_cn_plot(1:20:end), tau(1:16:end), u_cn_history(1:16:end, 1:20:end))
title('Crank-Nicolson')
xlabel('Asset price')
ylabel('Time to maturity')
zlabel('Call option value')
xlim([0, plot_limit])
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

function [x, u, tau, u_history] = crank_nicolson(K, r, sigma, T, h, Smax, Nt, plot_limit)
x = 0:h:Smax;
u = max(x - K, 0);

dt = T/Nt;
m = length(x);
plot_m = find(x <= plot_limit, 1, 'last');
tau = (0:Nt)*dt;
u_history = zeros(Nt + 1, plot_m);
u_history(1, :) = u(1:plot_m);

xi = x(2:end-1).';
diffusion = sigma^2*xi.^2/(2*h^2);
advection = r*xi/(2*h);
left = diffusion - advection;
center = -2*diffusion - r;
right = diffusion + advection;
N = length(xi);

A = sparse(1:N, 1:N, 1 - dt*center/2, N, N) ...
    + sparse(2:N, 1:N-1, -dt*left(2:N)/2, N, N) ...
    + sparse(1:N-1, 2:N, -dt*right(1:N-1)/2, N, N);
B = sparse(1:N, 1:N, 1 + dt*center/2, N, N) ...
    + sparse(2:N, 1:N-1, dt*left(2:N)/2, N, N) ...
    + sparse(1:N-1, 2:N, dt*right(1:N-1)/2, N, N);

for n = 1:Nt
    right_old = Smax - K*exp(-r*tau(n));
    right_new = Smax - K*exp(-r*tau(n + 1));

    rhs = B*u(2:end-1).';
    rhs(end) = rhs(end) + dt*right(end)*(right_old + right_new)/2;

    u(2:end-1) = (A\rhs).';
    u(1) = 0;
    u(end) = right_new;
    u_history(n + 1, :) = u(1:plot_m);
end
end

function price = black_scholes_call(S, tau, K, r, sigma)
d1 = (log(S/K) + (r + sigma^2/2)*tau)/(sigma*sqrt(tau));
d2 = d1 - sigma*sqrt(tau);
Phi = @(z) 0.5*erfc(-z/sqrt(2));
price = S.*Phi(d1) - K*exp(-r*tau).*Phi(d2);
end
