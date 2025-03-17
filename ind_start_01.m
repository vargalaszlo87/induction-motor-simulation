% Induction motor startup simulation (Frequency Ramp Only)
clear; clc; close all;

% Motor parameters
  V = 400; % Constant voltage [V]
  f_nom = 50; % Nominal frequency [Hz]
  p = 4; % Number of poles
  J = 9e-3; % Moment of inertia [kg·m^2]
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
  t_end = 1
  t = linspace(0, t_end, step); % 2 period
  id = zeros(size(t)); % D-tengely áram
  iq = zeros(size(t)); % Q-tengely áram
  idr = zeros(size(t)); % D-tengely rotor áram
  iqr = zeros(size(t)); % Q-tengely rotor áram
  s = zeros(size(t));
  tau = zeros(size(t)); % Nyomaték
  omega_r = zeros(size(t)); % Rotor szögsebessége  
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

  h = (t_end - 0) / step;
  for k = 2:length(t)
      % solving diff eq
      did_dt = (1 / Ls) * (Vsd(k) - Rs * id(k-1) + omega * Ls * iq(k-1));
      diq_dt = (1 / Ls) * (Vsq(k) - Rs * iq(k-1) - omega * Ls * id(k-1));
      
      % current refresh
      id(k) = id(k-1) + did_dt * h;
      iq(k) = iq(k-1) + diq_dt * h;
      
      % slip
      s(k) = max(0, min(1, (omega_s - (1 - s(k-1)) * omega_s) / omega_s));
      
      % Rotor differenciálegyenletek megoldása
      didr_dt = (1 / Lr) * (-Rr * idr(k-1) + s(k) * omega_s * Lm * iqr(k-1));
      diqr_dt = (1 / Lr) * (-Rr * iqr(k-1) - s(k) * omega_s * Lm * idr(k-1));

      % Rotor áramok frissítése
      idr(k) = idr(k-1) + didr_dt * h;
      iqr(k) = iqr(k-1) + diqr_dt * h;   
   
      % Elektromágneses nyomaték számítása
      tau(k) = (3/2) * p * Lm * (id(k) * iq(k) - idr(k) * iqr(k));

      % Rotor szögsebességének frissítése (helyes mozgásegyenlet!)
      omega_r(k) = omega_r(k-1) + (tau(k) / J) * h;   
      
      % n_rotor
      n_rotor(k) = (omega_r(k) / 2 * pi) * 60;
      
  end

  
  
% Ábrázolás
##figure;
##plot(t, Va, 'r', 'LineWidth', 1); hold on;
##plot(t, Vb, 'g', 'LineWidth', 1);
##plot(t, Vc, 'b', 'LineWidth', 1);
##xlabel('Idõ (s)'); ylabel('Feszültség (V)');
##title('Háromfázisú szinuszos feszültség');
##legend('V_a', 'V_b', 'V_c');
##grid on;;
