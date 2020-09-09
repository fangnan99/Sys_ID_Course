clc
disp('Running ps2 try.m');
clear

sampling_period_T = 1;
% number_cycles_r = 10;
samples_per_cycle_M = 1024;

j=1;
for number_cycles_r = [5,10,20,30,50]
    length_L = number_cycles_r * samples_per_cycle_M;
    gain_input = tf(0.1, [1, -1.7, 0.72], 1);
    gain_noise = tf([1, -0.92]*1.5, [1, -0.5], 1);

    input_u = [0; 2*repmat(randn(samples_per_cycle_M,1),number_cycles_r,1)];
    noise_e = 0.1*randn(length_L+1, 1);

    output_y = lsim(gain_input, input_u, 0:1:length_L) + lsim(gain_noise, noise_e, 0:1:length_L);

    % figure(1)
    % subplot(2,1,1)
    % plot(input_u)
    % subplot(2,1,2)
    % plot(output_y)

    autocorr_u = xcorr(input_u, length_L/2);
    % figure(2)
    % plot(-length_L/2:1:length_L/2, autocorr_u)

    average_period_u = input_u(2:1+samples_per_cycle_M);
    average_period_y = zeros(samples_per_cycle_M, 1);
    for i=2:1:number_cycles_r
        average_period_y = average_period_y + output_y((i-1)*samples_per_cycle_M+2:i*samples_per_cycle_M+1)/(number_cycles_r-1);
    end

    % figure(3)
    % subplot(2,1,1)
    % plot(average_period_u)
    % subplot(2,1,2)
    % plot(average_period_y)

    U = fft(average_period_u); % calculate N point FFTs
    Y = fft(average_period_y);
    omega = (2*pi/samples_per_cycle_M)*[0:samples_per_cycle_M-1]'; % frequency grid
    idx = find(omega > 0 & omega < pi); % positive frequencies

    % figure(4)
    % loglog(omega(idx),abs(U(idx)))
    % hold on
    % loglog(omega(idx),abs(Y(idx)))

    Gest = Y./U; % ETFE estimate
    Gfresp = squeeze(freqresp(gain_input,omega)); % "true" system response

    figure(j)
    subplot(2,1,1)
    loglog(omega(idx),abs(Gest(idx)))
    hold on
    loglog(omega(idx),abs(Gfresp(idx)))
    hold on
    subplot(2,1,2)
    semilogx(omega(idx),angle(Gest(idx)))
    hold on
    semilogx(omega(idx),angle(Gfresp(idx)))
    hold on
    Err = Gest - Gfresp; % calculate error
    subplot(2,1,1)
    hold on
    loglog(omega(idx),abs(Err(idx)))
    hold on
    
    j=j+1;
end
