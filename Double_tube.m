clc;
clear;
format short g

% -----------------------------
% Constants
% -----------------------------
pi = 3.141592653589793;
rhow = 999;            % Density of water (kg/m^3)
rho = 7800;            % Density of pipe material (kg/m^3)
g = 9.81;              % Gravity (m/s^2)
nut = 1.15e-6;         % Kinematic viscosity of water (m^2/s)
Q = 4.086;             % Flow rate (m^3/s)
e1 = 0.001;            % Base wall thickness (m)
k1 = 1.7;
k2 = 0.9;
e_toix = 0.2;          % Absolute roughness in mm
sep = 147099750;       % Max allowable stress (Pa)
L_hor = 3680;          % Horizontal length (m)
L_vert = 560;          % Vertical length (m)
h = 130;               % Hydraulic head (m)
en = 0.07;             % Energy selling price (€/kWh)

% -----------------------------
% Tables
% -----------------------------
htta = [0.819, 0.821, 0.829, 0.824, 0.832, 0.837, 0.840, ...
        0.842, 0.829, 0.83, 0.831, 0.831, 0.832, 0.832, ...
        0.832, 0.832, 0.832, 0.832] .* ...
       [0.97214, 0.9723, 0.97256, 0.97291, 0.97349, 0.97386, ...
       0.97406, 0.97416, 0.97416, 0.9742, 0.97425, 0.97428, ...
       0.97429, 0.9743, 0.97431, 0.97432, 0.97432, 0.97433];

cf = [0.5041, 0.4589, 0.4329, 0.4163, 0.397, 0.3868, 0.3807, ...
      0.3769, 0.3744, 0.3728, 0.3711, 0.3701, 0.3695, 0.369, ...
      0.3688, 0.3686, 0.3684, 0.3683];

D = [1066.8, 1117.6, 1168.4, 1219.2, 1320.8, ...
     1422.4, 1524.0, 1625.6, 1727.2, 1828.8, 1981.0, ...
     2134.0, 2286.0, 2438.0, 2591.0, 2743.0, 2896.0, 3048.0] * 1e-3;

emin = [7.14, 7.14, 7.14, 7.92, 7.92, 8.74, 8.74, 10.31, ...
        10.31, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, ...
        11.13, 11.13] * 1e-3;

% -----------------------------
% Preallocate arrays
% -----------------------------
n = length(D);
[es_hor, c_hor, Re_hor, l_hor, Dh_hor, ...
 Din, es_vert, c_vert, Re_vert, l_vert, Dh_vert, ...
 Pmax, e, ...
 CM_vert, CM_hor, CP_vert, CP_hor, CW_vert, CW_hor, ...
 f1_vert, f1_hor, f1, ...
 Nf_vert, Nf_hor, Ef_vert, Ef_hor, ...
 Cei_vert, Cei_hor, f2_vert, f2_hor, f2] = deal(zeros(1, n));

% -----------------------------
% Annuity factor (discounted sum)
% -----------------------------
annuityFactor = sum(1 ./ ((1 + 0.06).^(1:20)));

% -----------------------------
% Main loop
% -----------------------------
for k = 1:n
    % Horizontal pipe
    Din(k) = D(k) - 2 * emin(k); % Επειδή στην οριζόντια σωλήνωση h=0, το πάχος θα είναι το ελάχιστο δυνατόν
    es_hor(k) = e_toix / (Din(k) * 1000);
    c_hor(k) = (4 * Q) / (pi * Din(k)^2);
    Re_hor(k) = c_hor(k) * Din(k) / nut;

    % Friction factor (Colebrook)
    l_hor(k) = 1;
    for j = 1:50
        term = 2.51 / (Re_hor(k) * sqrt(l_hor(k))) + es_hor(k) / 3.71;
        l_new = (-2 * log10(term))^-2;
        if abs(l_new - l_hor(k)) < 1e-4
            break; 
        end
        l_hor(k) = l_new;
    end

    Dh_hor(k) = l_hor(k) * (L_hor / Din(k)) * (c_hor(k)^2 / (2 * g));
    Pmax(k) = 1.2 * (h - Dh_hor(k)) * rhow * g;

    % Wall thickness check
    e(k) = e1 + (D(k) * Pmax(k)) / (2 * (k2 / k1) * sep);
    if e(k) <= emin(k)
        e(k) = emin(k);
    end

    es_vert(k) = e_toix / (Din(k) * 1000);
    c_vert(k) = (4 * Q) / (pi * Din(k)^2);
    Re_vert(k) = c_vert(k) * Din(k) / nut;

    % Friction factor (Colebrook)
    l_vert(k) = 1;
    for j = 1:50
        term = 2.51 / (Re_vert(k) * sqrt(l_vert(k))) + es_vert(k) / 3.71;
        l_new = (-2 * log10(term))^-2;
        if abs(l_new - l_vert(k)) < 1e-4 
            break; 
        end
        l_vert(k) = l_new;
    end

    Dh_vert(k) = l_vert(k) * (L_vert / Din(k)) * (c_vert(k)^2 / (2 * g));

    % Cost calculation
    CM_vert(k) = rho * pi * D(k) * e(k) * L_vert * 1.16;
    CM_hor(k) = rho * pi * D(k) * emin(k) * L_hor * 1.16;
    CP_vert(k) = 32 * pi * D(k) * L_vert;
    CP_hor(k) = 32 * pi * D(k) * L_hor;
    CW_vert(k) = 1350 * D(k) * (L_vert / 6);
    CW_hor(k) = 1350 * D(k) * (L_hor / 6);

    f1_vert(k) = CM_vert(k) + CP_vert(k) + CW_vert(k);
    f1_hor(k) = CM_hor(k) + CP_hor(k) + CW_hor(k);
    f1(k) = f1_vert(k) + f1_hor(k);

    % Energy and economic return
    Nf_hor(k) = (rhow * g * Q * Dh_hor(k) * htta(k)) / 1000;
    Nf_vert(k) = (rhow * g * Q * Dh_vert(k) * htta(k)) / 1000;

    Ef_hor(k) = cf(k) * 8760 * Nf_hor(k);
    Ef_vert(k) = cf(k) * 8760 * Nf_vert(k);

    Cei_hor(k) = en * Ef_hor(k);
    Cei_vert(k) = en * Ef_vert(k);

    f2_hor(k) = Cei_hor(k) * annuityFactor;
    f2_vert(k) = Cei_vert(k) * annuityFactor;
    f2(k) = f2_vert(k) + f2_hor(k);
end

% -----------------------------
% Plot results
% -----------------------------
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

pipe_indices = 23:40;  % Standard pipe numbering

plot(pipe_indices, (f1_vert + f2_vert) / 1e6, '-', 'LineWidth', 2, ...
     'Marker', 'o', 'MarkerSize', 4, ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
     'Color', colors);

title('Συνολικό κόστος σωλήνωσης και απωλειών συναρτήσει της διαμέτρου για το κάθετο τμήμα ');
xlabel('Αριθμός τυποποίησης σωλήνωσης [-]');
ylabel('Συνολικό κόστος κάθετη σωλήνωσης [εκατομμύρια €]');
xlim([23 40])

% -----------------------------
% Optimal diameter
% -----------------------------
[~, index] = min(f1 + f2);
fprintf('Βέλτιστη διάμετρος: %.2f mm\n', D(index) * 1000);
k_hor = (8 * l_hor(index) * L_hor) / (pi^2 * D(index)^5 * g);
k_vert = (8 * l_vert(index) * L_vert) / (pi^2 * D(index)^5 * g);
