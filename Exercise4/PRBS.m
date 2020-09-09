clear
clc

figure(1)
for order = 5:1:8
    signal = idinput(2^order-1, 'prbs', [0,1], [-5,5]);
    length = size(signal,1);
    omega = (2*pi/length) * [0:length-1]';
    idx = find(omega >= 0 & omega <= pi);
    spectral = mypdg(fft(signal));
    loglog(omega(idx), spectral)
    hold on
end