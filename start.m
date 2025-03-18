close all;
clear;
clc;

global w_rotor = 0;

% LOCALE FUNCTION

  function dXdt = motor_dynamics(t, X, params)


    % expand parameters

      Rr = params(1);
      R_double_prime = params(2);
      sigma_Ls = params(3);
      Lr = params(4);
      Lm = params(5);
      w_s = params(6);  % Szinkron sebesség (hálózati frekvencia)
      w_r = params(7);
      V = params(8);
      J = params(9);  % Rotor tehetetlenségi nyomaték
      B_ = params(10);  % Súrlódási együttható
      P = params(11); % Póluspárok száma
      T_L = params(12); % Terhelő nyomaték

    % input signal

      % three pahse signal
      ua = V * cos(w_s * t);
      ub = V * cos(w_s * t - 2*pi/3);
      uc = V * cos(w_s * t + 2*pi/3);

      % Clarke transformation
      u_alpha = (2/3) * (ua - 0.5*ub - 0.5*uc);
      u_beta  = (2/sqrt(3)) * (ub - uc);

      % Park transformation
      u_sd = u_alpha * cos(w_s * t) + u_beta * sin(w_s * t);
      u_sq = -u_alpha * sin(w_s * t) + u_beta * cos(w_s * t);

    % expand state vector
      psi_rd = X(1);
      i_sd  = X(2);
      i_sq  = X(3);

    % matrices of differential equations

      A = [
          -Rr / Lr, (Lm * Rr) / Lr, 0;
          (Lm * Lr) / (pow2(Lr) * sigma_Ls), -R_double_prime / sigma_Ls, w_s;
          -w_r * Lm / (Lr * sigma_Ls), -w_s, -R_double_prime / sigma_Ls
      ];


      B = [
          0, 0;
          1 / sigma_Ls, 0;
          0, 1 / sigma_Ls
      ];

      U = [
          u_sd;
          u_sq
      ];

    % calc
    dXdt = [A * X + B * U];

  end;


% MAIN SECTION

  % motor parameters

    Rr = 0.5;
    Lr = 0.70315;
    Lm = 0.19779;
    Ls =  0.7205;
    Rs = 1.5396;
    P = 4;  % polusparok szama
    J = 0.01;  % Tehetetlensegi nyomatek
    B = 0.002;  % Surloadsei egyutthato
    T_L = -2;  % Terhelo nyomatek (Nm)

  % simulation parameters

    % signal
    A = 400;
    f = 50;
    w_s = 2 * pi * f;

    % time
    t_start = 0;
    t_end = 2;
    step = 1e4;
    tspan = linspace(t_start, t_end, step);

  % calculated parameters

    V = A * sqrt(2);

  % differential equations parameteres

    % auxiliary function / parameter
    w_r = 0;
    sigma_Ls = Ls - (Lm^2 / Lr);
    R_double_prime = (Lm^2 / Lr^2) * Rr + Rs;


% Params vector

    % params vector

    %params = [Rr, R_double_prime, sigma_Ls, Lr, Lm, w_rotor, w, V]
    params = [Rr, R_double_prime, sigma_Ls, Lr, Lm, w_s, w_r, V, J, B, P, T_L];

  % initial condition: psi_rd, i_sd, i_sq

    X0 = [0; 0; 0];

  % ode solving

    [t, X] = ode45(@(t, X) motor_dynamics(t, X, params), tspan, X0);


% OUTPUT SECTION

  figure;

  subplot(5,1,1);
  plot(t, X(:,1));
  xlabel('Time (s)'); ylabel('F RD');

  subplot(5,1,2);
  plot(t, X(:,2));
  xlabel('Time (s)'); ylabel('I SD');

  subplot(5,1,3);
  plot(t, X(:,3));
  xlabel('Time (s)');
  ylabel('I SQ');


  grid on;










