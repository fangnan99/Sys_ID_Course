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

% Generate and plot time domain windows
figure(4)
hold on
for gamma = [5, 10, 50, 100]
    [lags,WHann] = WtHann(gamma,N);
    plot(lags, WHann) ; 
end
legend ('5', '10', '50', '100') ;

figure (5) ;
loglog(omega, abs(Gfresp))
for gamma = [5, 10, 50, 100]
    j = 0 ;
    [lags,WHann] = WtHann(gamma,N);
    autocorrelation_u = 0*WHann;
    crosscorrelation_u = 0*WHann;
    for tao = -(N/2 -1):(N/2)
        j = j+1 ;
        autocorrelation = 0;
        for k = 1:size(lags,1)
            ktao = k-tao;
            while (ktao<=0) 
                ktao = ktao + N;
            end
            while (ktao > N)
                ktao = ktao - N;
            end
            autocorrelation = autocorrelation + input_u(k)*input_u(ktao);
        end
        autocorrelation_u(j) = autocorrelation;
    end
    denominator = fft( WHann .* autocorrelation_u);
    
    j = 0;
    for tao = -(N/2 -1):(N/2)
        j = j+1;
        crosscorrelation = 0;
        for k = 1:N
            ktao = k-tao;
            while (ktao<=0) 
                ktao = ktao + N;
            end
            while (ktao > N)
                ktao = ktao - N;
            end
            crosscorrelation = crosscorrelation + output_y(k)*input_u(ktao);
        end
        crosscorrelation_u(j) = crosscorrelation;
    end
    numerator = fft( WHann .* crosscorrelation_u); 
    
    %compute smooth transfer function
    G_smooth_time = abs(numerator./denominator);

    %plot results
    hold on
    loglog(omega, abs(G_smooth_time));
    hold on;
    title('time domain window')
end 
legend('True', '5', '10', '50', '100')