% Solution to Problem Set 5 MATLAB Exercise
clear

% Provided Parameters
zeta_z = 0.1;
omega_z = 3.0;
zeta_p = 0.1;
omega_p = 3.5;
T_s = 0.02;
noise_var = 0.1;
% Defined Parameters


% Define transfer functions and convert to discrete time domain
G_s = tf(conv([1 2*zeta_z*omega_z omega_z^2], [5000]), conv([1 2*zeta_p*omega_p omega_p^2], conv([1 50], [1 200])));
G_dz = c2d(G_s, T_s, 'zoh');
C_dz = tf([1.25 -0.75], [1 -1], T_s);

%% 1. 
% Calculate S_dz

G_dz_nu = poly2sym(G_dz.Numerator{1});
G_dz_de = poly2sym(G_dz.Denominator{1});
C_dz_nu = poly2sym(C_dz.Numerator{1});
C_dz_de = poly2sym(C_dz.Denominator{1});

S_dz_de = poly2sym([1]) + G_dz_nu * C_dz_nu / (G_dz_de * C_dz_de);
S_dz_sym = poly2sym([1]) / S_dz_de;

S_dz_sym = simplify(S_dz_sym);
[S_dz_nu S_dz_de] = numden(S_dz_sym);
S_dz = tf(sym2poly(S_dz_nu), sym2poly(S_dz_de), T_s);

% S_dz = 1/(1-G_dz*C_dz);

% Calculate and show poles of S_dz
disp('The poles of S_dz are');
disp(pole(S_dz));
disp('Their absolute values are');
disp(abs(pole(S_dz)));

% Calculate T_dz = 1- S_dz
T_dz_sym = 1-S_dz_sym;
[T_dz_nu T_dz_de] = numden(T_dz_sym);
T_dz = tf(sym2poly(T_dz_nu), sym2poly(T_dz_de), T_s);

%% 2.
N_spp = 1023;    % No samples per period
N_p = 6;        % No Periods
r_period = idinput(N_spp, 'prbs', [0 1], [-1, 1]);
r = repmat(r_period, N_p, 1);
v = 1*randn(N_spp*N_p, 1);
y = lsim(T_dz, r) + lsim(S_dz, v);
u = lsim(C_dz*S_dz, r-v);
figure(1)
plot(u);
hold on
plot(y);

omega = (2*pi/N_spp)*[0:N_spp-1]'; % frequency grid
omega_z = exp(i*omega);
idx = find(omega > 0 & omega < pi); % positive frequencies

Gfresp = squeeze(freqresp(G_dz,omega_z)); % True system frequency response
y = reshape(y, N_spp, N_p);
y_period = mean(y(:,2:N_p), 2);
u = reshape(u, N_spp, N_p);
u_period = mean(u(:,2:N_p), 2);

U_period = fft(u_period);
Y_period = fft(y_period);

Getfe = Y_period./U_period;

figure(2)
loglog(omega(idx)*50, abs(Gfresp(idx)))
hold on
loglog(omega(idx)*50, abs(Getfe(idx)))
hold on

%% 3.
e = r - reshape(y, N_spp*N_p, 1);
e = reshape(e, N_spp, N_p);
e_period = mean(e(:,N_p), 2);

R_period = fft(r_period);
E_period = fft(e_period);

Setfe = E_period./R_period;
Sfresp = squeeze(freqresp(S_dz,omega_z));

figure(3)
loglog(omega(idx)*50, abs(Sfresp(idx)));
hold on
loglog(omega(idx)*50, abs(Setfe(idx)));

%% 4.
w_period = idinput(N_spp, 'prbs', [0 1], [-1, 1]);
w = repmat(w_period, N_p, 1);
y = lsim(G_dz*S_dz, w) + lsim(S_dz, v);
u = w - lsim(C_dz*S_dz, v);

GSfresp = squeeze(freqresp(G_dz*S_dz,omega_z));

w = reshape(w, N_spp, N_p);
y = reshape(y, N_spp, N_p);
y_period = mean(y(:,2:N_p),2);

W_period = fft(w_period);
Y_period = fft(y_period);

GSetfe = Y_period./W_period;

figure(4)
loglog(omega(idx)*50, abs(GSfresp(idx)));
hold on
loglog(omega(idx)*50, abs(GSetfe(idx)));

%% 5.
figure(2)
hold on
loglog(omega(idx)*50, abs(GSetfe(idx)./Setfe(idx)))