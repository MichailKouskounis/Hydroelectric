clc;
clear;
format short g

% Constants
pi = 3.14159265358979323846;
rhow = 999;         % kg/m^3
rho = 7800;         % kg/m^3 (steel)
g = 9.81;           % m/s^2
nut = 1.15e-6;      % m^2/s (kinematic viscosity)
Q = 4.086;          % m^3/s

% Efficiency and cost factor data (length 19)
htta = [0.818, 0.821, 0.828, 0.824, 0.828, 0.835, 0.839, 0.841, 0.843, 0.83, ...
        0.83, 0.831, 0.831, 0.832, 0.832, 0.832, 0.832, 0.832, 0.832] .* ...
       [0.9721, 0.97228, 0.97251, 0.97289, 0.97323, 0.9737, 0.974, 0.97412, ...
        0.97421, 0.97419, 0.97422, 0.97426, 0.97429, 0.9743, 0.97431, ...
        0.97432, 0.97432, 0.97432, 0.97433];

cf = [0.5150, 0.4636, 0.435, 0.4173, 0.4054, 0.3911, 0.3831, 0.3784, 0.3753, ...
       0.3733, 0.372, 0.3706, 0.3698, 0.3692, 0.3689, 0.3686, 0.3685, ...
       0.3684, 0.3683];

% Pipe outer diameters (m) and minimum wall thicknesses (m)
D = [1016.0, 1066.8, 1117.6, 1168.4, 1219.2, 1320.8, ...
     1422.4, 1524.0, 1625.6, 1727.2, 1828.8, 1981.0, ...
     2134.0, 2286.0, 2438.0, 2591.0, 2743.0, 2896.0, 3048.0] * 1e-3;

emin = [7.14, 7.14, 7.14, 7.14, 7.92, 7.92, 8.74, 8.74, ...
        10.31, 10.31, 11.13, 11.13, 11.13, 11.13, 11.13, ...
        11.13, 11.13, 11.13, 11.13] * 1e-3;

% Design parameters
e1 = 0.001;         % m
k1 = 1.7;
k2 = 0.9;
e_toix = 0.2;       % mm
sep = 147099750;    % Pa
L = 3870;           % m
h = 130;            % m
en = 0.07;          % €/kWh

% Max pressure
Pmax = (1.2 * h) * rhow * g;

% Initialization
n = length(D);
f1 = zeros(1, n);
f2 = zeros(1, n);
e = zeros(1, n);
Din = zeros(1, n);
Nf = zeros(1, n);
Ef = zeros(1, n);
l = zeros(1, n);
Dh = zeros(1,n);
CM = zeros(1,n);
CP = zeros(1,n);
CW = zeros(1,n);

% Discount factor sum for 20 years at 6% interest
discount_sum = sum(1 ./ ((1 + 0.06).^(1:20)));

% Loop over each diameter
for k = 1:n
    % Wall thickness calculation
    e(k) = e1 + (D(k) * Pmax) / (2 * (k2 / k1) * sep);
    if e(k) <= emin(k)
        e(k) = emin(k);
    end

    Din(k) = D(k) - 2 * e(k);

    % Cost function f1
    CM(k) = rho * pi * D(k) * e(k) * L * 1.16;
    CP(k) = 32 * pi * D(k) * L;
    CW(k) = 1350 * D(k) * (L / 6);
    f1(k) = CM(k) + CP(k) + CW(k);

    % Efficiency-related calculations
    es = e_toix / (Din(k) * 1000);  % relative roughness
    c = (4 * Q) / (pi * Din(k)^2); % velocity
    Re = (c * Din(k)) / nut;

    % Colebrook equation
    l_old = 0.02;
    l(k) = 0.03;
    for j = 1:50
        term = 2.51 / (Re * sqrt(l(k))) + es / 3.71;
        l_new = (-2 * log10(term))^-2;
        if abs(l_new - l_old) < 1e-4
            break;
        end
        l_old = l_new;
        l(k) = l_new;
    end

    % Head loss and energy recovery
    Dh(k) = l(k) * (L / Din(k)) * (c^2 / (2 * g));
    Nf(k) = (rhow * g * Q * Dh(k) * htta(k)) / 1000; % kW
    Ef(k) = cf(k) * 8760 * Nf(k);                % kWh/year
    Cei = en * Ef(k);                            % €/year

    f2(k) = Cei * discount_sum;
end

% Plotting
figure;
hold on
set(gca,'Yscale','log')
set(gcf, 'PaperPositionMode', 'auto'); 
grid on;
grid minor;
ax.GridColor = [0.85 0.85 0.85];
ax.GridAlpha = 0.7;
ax.MinorGridColor = [0.9 0.9 0.9];
ax.MinorGridLineStyle = ':';
set(gca, 'FontSize', 14);
%set(gca, 'YScale', 'log');  % Not log
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 5]);
ax=gca;
colors = [0.5 0.3 0.5]; % Color for plot
% Plot total cost
pipe_indices = 22:40;  % Standard pipe numbering
plot(pipe_indices, (f1 + f2)/1e6, '-', 'LineWidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
    'Color', [0.5 0.3 0.5]);

title('Συνολικό κόστος σωλήνωσης και απωλειών συναρτήσει της διαμέτρου');
xlabel('Αριθμός τυποποίησης σωλήνωσης [-]');
ylabel('Συνολικό κόστος σωλήνωσης [εκατομμύρια €]');

% Display optimal diameter
[~, index] = min(f1 + f2);
optimal_D = D(index);

fprintf('Βέλτιστη διάμετρος D = %.3f m\n', optimal_D);

% Compute and display coefficient k
k_coef = (8 * l(index) * L) / (pi^2 * Din(index)^5 * g);
fprintf('k = %.3e\n', k_coef);
