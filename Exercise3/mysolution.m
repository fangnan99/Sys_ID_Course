clear
clc

% Define transfer functions
gain_input = tf([0.1, 0], conv([1, -1.7, 0.72], [1, -0.98, 0.9]), 1);
gain_noise = tf([1, -0.9]*0.5, [1, -0.25], 1);

% Generate signal
num_samples = 1024;
N = 1024;
input_u = randn(num_samples, 1);
input_noise_e = randn(num_samples, 1);

% Run simulation
output_y = lsim(gain_input, input_u, 0:1:num_samples-1) + lsim(gain_noise, input_noise_e, 0:1:num_samples-1);

% Plot signal vs time
figure(1)
subplot(2,1,1)
plot(input_u)
hold on
plot(input_noise_e)
subplot(2,1,2)
plot(output_y)

% Calculate FFT
U = fft(input_u); % calculate N point FFTs
Y = fft(output_y);

% Frequency grid for full length estimations
omega = (2*pi/num_samples)*[0:num_samples-1]'; % frequency grid
idx = find(omega > 0 & omega < pi); % positive frequencies

% Plot FFT
figure(2)
loglog(omega(idx),abs(U(idx)))
hold on
loglog(omega(idx),abs(Y(idx)))

% ETFE & True system response
Gest = Y./U; % ETFE estimate
Gfresp = squeeze(freqresp(gain_input,omega)); % "true" system response
Hfresp = squeeze(freqresp(gain_noise,omega));

% Plot ETFE & True system response
figure(3)
subplot(2,1,1)
loglog(omega(idx),abs(Gest(idx)))
hold on
loglog(omega(idx),abs(Gfresp(idx)))
subplot(2,1,2)
semilogx(omega(idx),angle(Gest(idx)))
hold on
semilogx(omega(idx),angle(Gfresp(idx)))
Err = Gest - Gfresp; % calculate error
subplot(2,1,1)
% hold on
% loglog(omega(idx),abs(Err(idx)))
hold on
loglog(omega(idx),abs(Hfresp(idx)))

% Split data into 4 chunks and perform estimation on each chunk
n_split = 4;
Gest_all = zeros(num_samples/n_split, n_split);
for i=1:1:n_split
    part_u = input_u((i-1)*num_samples/n_split+1:i*num_samples/n_split);
    part_y = output_y((i-1)*num_samples/n_split+1:i*num_samples/n_split);
    part_U = fft(part_u);
    part_Y = fft(part_y);
    Gest_all(:,i) = part_Y./part_U;
end
Gest_avg = mean(Gest_all, 2);

% Frequency grid for chunked estimation
omega_avg = (2*pi/(num_samples/n_split))*[0:num_samples/n_split-1]'; % frequency grid
idx_avg = find(omega_avg > 0 & omega_avg < pi); % positive frequencies

% Plot
figure(4)
subplot(2,1,1)
loglog(omega_avg(idx_avg),abs(Gest_avg(idx_avg)))
hold on
loglog(omega(idx),abs(Gfresp(idx)))
subplot(2,1,2)
semilogx(omega_avg(idx_avg),angle(Gest_avg(idx_avg)))
hold on
semilogx(omega(idx),angle(Gfresp(idx)))
% Err = Gest_avg - Gfresp; % calculate error
subplot(2,1,1)
% hold on
% loglog(omega(idx),abs(Err(idx)))
% hold on
% loglog(omega(idx),abs(Hfresp(idx)))
hold on
loglog(omega(idx),abs(Gest(idx)))

% Generate and plot f domain Hann Windows
figure(5)
for gamma = [5, 10, 50, 100]
    [omega, hann_win_f] = WfHann(gamma, N);
    hold on
    plot(hann_win_f)
    hold on
end

% Windowed estimation, code from lecture notes
Gs = 0*Gest;
figure(6)
loglog(omega(idx),abs(Gfresp(idx)))
hold on
for gamma = [5, 10, 50, 80]
    [omega, Wg] = WfHann(gamma, N);
    zidx = find(omega == 0);
    omega = [omega(zidx:N); omega(1:zidx-1)];
    Wg = [Wg(zidx:N), Wg(1:zidx-1)];
    a = U.*conj(U);
    for wn = 1:N
        Wnorm = 0;
        for xi = 1:N
            widx = mod(xi-wn,N)+1;
            Gs(wn) = Gs(wn) + Wg(widx) * Gest(xi) * a(xi);
            Wnorm = Wnorm + Wg(widx) * a(xi);
        end
        Gs(wn) = Gs(wn)/Wnorm;
    end
    loglog(omega, abs(Gs))
    hold on
end
legend('True', 'Gamma=5', 'Gamma=10', 'Gamma=50', 'Gamma=80')

