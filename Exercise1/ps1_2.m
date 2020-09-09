% Problem here:

% What is the relationship between periodogram furmula in problemset &
% matlab periodogram function.
% Observed difference: Matlab result is 3.14 * calculated result, while
% 6.28* for 1st and last data.

% What is the defination of the gain obtained from bode().
% why 20*log10(bode_magnitude) = 10*log10(avg_calculated result)

plant = tf(1, [1, -0.9, 0.5], 1);
total_y = zeros(513, 1);

round = 5;
for a=1:1:round
    t = 1024;
    e_k = randn(t, 1);
    E_w = mydft(e_k);
    w = 0:2*pi/1024:pi;
    e_pxx_ml = periodogram(e_k);
    e_my_pxx = mypdg(E_w);
   
    y_k = lsim(plant, e_k, 0:1:t-1);

    [bodemag, bodephase, wout]=bode(plant, w);

    Y_w = mydft(y_k);
    y_my_pxx = mypdg(Y_w);
    total_y = total_y + y_my_pxx;
    disp(a);
end
avg_y = total_y/round;
mag = reshape(bodemag, 513, 1);
plot(log10(mag)*20);
hold on;
plot(10*log10(avg_y));