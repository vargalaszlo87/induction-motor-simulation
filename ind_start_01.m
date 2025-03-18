% Induction motor startup simulation (Frequency Ramp Only)

clear; clc; close all;

% Motor parameters

  V = 400; % Constant voltage [V]
  f_nom = 50; % Nominal frequency [Hz]
  p = 4; % Number of poles
  J = 9e-3; % Moment of inertia [kg�m^2]
  M = 4.89e-2 % Mutual inductance [H]
  Rs = 1.5293; % Stator resistance [?]
  Rr = 0.7309; % Rotor resistance [?]
  Lm = 0.19779; % Magnetizing inductance [H]
  Lls = 3.56e-3; % Stator leakage inductance [H]
  Llr = 5.35e-3; % Rotor leakage inductance [H]
  Lr = 0.70315; % Rotor inductance [H]
  Ls = 0.72323; % Stator inductance [H]
  omega_s = 2 * pi * f_nom; % sync omega [rad/s]

% 50Hz

  f = 50
  A = 400 * sqrt(2);
  omega = 2 * pi * f;

% Time simulation setting

  step = 10e3;
  %t_end = 1/f * 50;
  t_end = 1;
  t = linspace(0, t_end, step); % 2 period
  dt = t_end / step;
  id = zeros(size(t)); % D-tengely �ram
  iq = zeros(size(t)); % Q-tengely �ram
  idr = zeros(size(t)); % D-tengely rotor �ram
  iqr = zeros(size(t)); % Q-tengely rotor �ram
  s = zeros(size(t));
  tau = zeros(size(t)); % Nyomat�k
  omega_r = zeros(size(t)); % Rotor sz�gsebess�ge
  n_rotor = zeros(size(t));

% Three phase (stator)

  Vsa = A * sin(omega * t);
  Vsb = A * sin(omega * t - 2 * pi / 3);
  Vsc = A * sin(omega * t + 2 * pi / 3);

% convert a,b,c to alfa, beta (Clarke)

  Vs_alpha = (2/3) * (Vsa - 0.5 * Vsb - 0.5 * Vsc);
  Vs_beta  = (2/sqrt(3)) * (0.5 * (Vsb - Vsc));

% convert alfa, beta to d, q (Park)

  theta = omega * t;
  Vsd = Vs_alpha .* cos(theta) + Vs_beta .* sin(theta);
  Vsq = -Vs_alpha .* sin(theta) + Vs_beta .* cos(theta);

% stator current (Euler)


% Differenci�legyenletek Runge-Kutta megold�ssal
for k = 2:length(t)
    % Cs�sz�s sz�m�t�sa
    s(k) = max(0, min(1, (omega_s - omega_r(k-1)) / omega_s));

    % Stator differenci�legyenletek (Runge-Kutta 4. rend� m�dszer)
    k1_id = (1 / Ls) * (Vsd(k) - Rs * id(k-1) + omega_s * Ls * iq(k-1));
    k1_iq = (1 / Ls) * (Vsq(k) - Rs * iq(k-1) - omega_s * Ls * id(k-1));

    k2_id = (1 / Ls) * (Vsd(k) - Rs * (id(k-1) + k1_id*dt/2) + omega_s * Ls * (iq(k-1) + k1_iq*dt/2));
    k2_iq = (1 / Ls) * (Vsq(k) - Rs * (iq(k-1) + k1_iq*dt/2) - omega_s * Ls * (id(k-1) + k1_id*dt/2));

    id(k) = id(k-1) + k2_id * dt;
    iq(k) = iq(k-1) + k2_iq * dt;

    % Rotor differenci�legyenletek (Runge-Kutta 4)
    k1_idr = (1 / Lr) * (-Rr * idr(k-1) + s(k) * omega_s * Lm * iqr(k-1));
    k1_iqr = (1 / Lr) * (-Rr * iqr(k-1) - s(k) * omega_s * Lm * idr(k-1));

    k2_idr = (1 / Lr) * (-Rr * (idr(k-1) + k1_idr*dt/2) + s(k) * omega_s * Lm * (iqr(k-1) + k1_iqr*dt/2));
    k2_iqr = (1 / Lr) * (-Rr * (iqr(k-1) + k1_iqr*dt/2) - s(k) * omega_s * Lm * (idr(k-1) + k1_idr*dt/2));

    idr(k) = idr(k-1) + k2_idr * dt;
    iqr(k) = iqr(k-1) + k2_iqr * dt;

    % Elektrom�gneses nyomat�k sz�m�t�sa
    tau(k) = (3/2) * p * Lm * (id(k) * iq(k) - idr(k) * iqr(k));

    % Rotor sz�gsebess�g�nek friss�t�se
    omega_r(k) = omega_r(k-1) + (tau(k) / J) * dt;
end

% �br�zol�s
figure;
subplot(4,1,1);
plot(t, id, 'r', t, iq, 'b', 'LineWidth', 1.5);
xlabel('Id� (s)'); ylabel('Stator �ram (A)');
legend('i_d', 'i_q');
title('Stator �ramok');

subplot(4,1,2);
plot(t, idr, 'r', t, iqr, 'b', 'LineWidth', 1.5);
xlabel('Id� (s)'); ylabel('Rotor �ram (A)');
legend('i_{dr}', 'i_{qr}');
title('Rotor �ramok');

subplot(4,1,3);
plot(t, s, 'k', 'LineWidth', 2);
xlabel('Id� (s)'); ylabel('Cs�sz�s');
title('Cs�sz�s id�f�ggv�nye');

subplot(4,1,4);
plot(t, tau, 'm', 'LineWidth', 2);
xlabel('Id� (s)'); ylabel('Nyomat�k (Nm)');
title('Elektrom�gneses nyomat�k');

grid on;



