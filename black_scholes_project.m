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

%% Article method and simple verification
exact_price = black_scholes_call(S0, T, K, r, sigma);
paper_price = [3.1676, 3.7501, 3.8871, 3.9189];
price = zeros(size(Nt_values));
error = zeros(size(Nt_values));

for q = 1:length(Nt_values)
    Nt = Nt_values(q);
    [x, u] = ade_method(K, r, sigma, T, h, Smax, Nt, -r);

    [~, target_index] = min(abs(x - S0));
    price(q) = u(target_index);
    error(q) = abs(price(q) - exact_price);
end

observed_order = [NaN, log2(error(1:end-1)./error(2:end))];
difference_from_paper = abs(price - paper_price);

results = table(Nt_values', price', paper_price', exact_price*ones(length(Nt_values),1), ...
    error', observed_order', difference_from_paper', ...
    'VariableNames', {'Nt','price','paper_price','exact_price','error', ...
    'observed_order','difference_from_paper'});
disp(results)

% The paper prints four decimals, so half a last printed digit is the target.
assert(max(difference_from_paper) < 5e-4, 'The paper prices were not reproduced.')
assert(all(diff(error) < 0), 'The error did not decrease under time refinement.')
assert(error(end) < 2e-2, 'The finest result is not close to the exact price.')
assert(abs(u(1)) < 1e-12 && all(isfinite(u)), 'Boundary or finite-value check failed.')

% A wrong sign on the reaction term must be detected by the exact solution.
[x_wrong, u_wrong] = ade_method(K, r, sigma, T, h, Smax, Nt_values(end), r);
[~, target_index] = min(abs(x_wrong - S0));
wrong_error = abs(u_wrong(target_index) - exact_price);
assert(wrong_error > error(end) && abs(u_wrong(target_index) - price(end)) > 2e-2, ...
    'The verification did not detect the planted sign error.')

fprintf('All simple checks passed.\n')

%% Plot the finest numerical solution
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
function [x, u] = ade_method(K, r, sigma, T, h, Smax, Nt, reaction)
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
            + reaction*u(i)/2;
        denominator = 1/dt + diffusion - advection - reaction/2;
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
            + reaction*u(i)/2;
        denominator = 1/dt + diffusion + advection - reaction/2;
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
