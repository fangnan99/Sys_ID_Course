function [p1_R,p1_omega,p1_a,p1_var] = HS2019_SysID_midterm_p1_11235813()

%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);
LegiNumber='11235813';
p1_U = HS2019_SysID_midterm_p1_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variable p1_U to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Task 1: Calculation of Autocorrelation

fig = 1;    % figure number
[S,K] = size(p1_U);
Ts = 0.5;

% Identify period using two methods
% method 1: plot inputs in time domain and identify the period T
figure(fig); fig=fig+1;
for s = 1:S
    subplot(S,1,s)
    plot(p1_U(s,:))
    hold on 
    title(['u signal, sensor ', num2str(s)])
    ylim([min(p1_U(1,:))-10, max(p1_U(1,:))+10])
end

% method 2: compute autocorrelation of input to identify the period T    
L =K;
tau = -L/2+1:L/2;
Ru = NaN*ones(S,length(tau));
for s = 1:S
    for t = 1:length(tau)
        Ru(s,t) = 0;
        for k = 1:L
            idx = k-tau(t);
            if idx <1
                idx = k-tau(t)+L;
            elseif idx >L
                idx = k-tau(t)-L;
            end
            Ru(s,t) = Ru(s,t) + 1/(L) *p1_U(s,k)*p1_U(s,idx);       
        end
    end
end
% Compute mean to average out the noise
Ru_mean = mean(Ru(:,:),1);

figure(fig); fig=fig+1;
plot(tau,Ru_mean)
xlabel('\tau')
title('Autocorrelation of the whole signal')
ylabel('R_u(\tau)')

% Identify true period
T = 240;  % by inspecting the autocorrelation plot (or the time domain measurements)
N = 6; % Identify number of periods in the signal from plots

% Compute the autocorrelation of first row, first period using the period length T
L= T;
tau = -L/2+1:L/2;
p1_R = NaN*ones(1,length(tau));
for t = 1:length(tau)
    p1_R(1,t) = 0;
    for k = 1:L
        idx = k-tau(t);
        if idx <1
            idx = k-tau(t)+L;
        elseif idx >L
            idx = k-tau(t)-L;
        end
        p1_R(1,t) = p1_R(1,t) + 1/(L) *p1_U(1,k)*p1_U(1,idx);      
    end
end

% Plot autocorrelation for one period
figure(fig); fig=fig+1;
plot(tau,p1_R)
xlabel('\tau')
title('Autocorrelation of one period')
ylabel('R_u(\tau)')
%% Task 2: Estimation of signal components

U_dft = NaN*ones(S,N*T);
Psd_U = NaN*ones(S,N*T);

omegas = 2*pi/(N*T)/Ts*[0:N*T-1];

% Compute average spectrum
for s = 1:S
    U_dft(s,:) = fft(p1_U(s,1:N*T));
    Psd_U(s,:) = 1/(N*T) * abs(U_dft(s,1:N*T)).^2;
end
Psd_u = mean(Psd_U,1);

% Plot average spectrum to identify peaks
figure(fig); fig=fig+1;
plot(omegas,Psd_u, '*')
xlabel('\omega [rad/sec]')
title('Average PSD')

% Identify right level to cut off signals (100-280)
idx_greater = find(Psd_u > 100 );
omegas_excited = omegas(idx_greater);  % excited frequencies in [rad/sec]  (in increasing order)

p1_omega = omegas_excited(1:length(omegas_excited)/2);  % only frequencies before pi/T
M = length(p1_omega); 

input_energies = Psd_u(idx_greater);
input_energies = input_energies(1:M);
p1_a = NaN*ones(1,M);
for i = 1:M
    % compute alpha(i) coefficient    
    % each sinusoid has energy =  N*T * alpha(i)^2/2 = input_frequencies(i)*2 (because N*T has an integer num of
    % periods of the sinusoid) 
    p1_a(i) = sqrt(2/(N*T)*    input_energies(i)*2);
end


%% Task 3: Estimation of noise variance

% Identify frequencies not in the underlying signal and find the mean of
% the Psd at these frequencies. The noise has a flat spectrum

p1_var = mean(Psd_u(~ismember(omegas,omegas_excited)));

end

