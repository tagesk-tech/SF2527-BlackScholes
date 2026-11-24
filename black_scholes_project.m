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

%% Article method
[x, u] = ade_method(K, r, sigma, T, h, Smax, Nt);

[~, target_index] = min(abs(x - S0));
numerical_price = u(target_index);
exact_price = black_scholes_call(S0, T, K, r, sigma);
fprintf('Numerical price: %.6f\n', numerical_price)
fprintf('Exact price:     %.6f\n', exact_price)

%% Compare with the exact solution
exact_profile = black_scholes_call(x, T, K, r, sigma);
show = x <= 160;

figure;
plot(x(show), exact_profile(show), 'k-', 'LineWidth', 1.5)
hold on
plot(x(show), u(show), '--', 'LineWidth', 1.3)
hold off
xlabel('Asset price')
ylabel('Call option value')
legend('Exact', 'ADE method', 'Location', 'northwest')
grid on

%% Functions
function [x, u] = ade_method(K, r, sigma, T, h, Smax, Nt)
x = 0:h:Smax;
u = max(x - K, 0);
phi = zeros(size(u));
psi = zeros(size(u));

dt = T/Nt;
m = length(x);

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
    m = mnew;
end

x = x(1:m);
u = u(1:m);
end

function price = black_scholes_call(S, tau, K, r, sigma)
d1 = (log(S/K) + (r + sigma^2/2)*tau)/(sigma*sqrt(tau));
d2 = d1 - sigma*sqrt(tau);
Phi = @(z) 0.5*erfc(-z/sqrt(2));
price = S.*Phi(d1) - K*exp(-r*tau).*Phi(d2);
end
