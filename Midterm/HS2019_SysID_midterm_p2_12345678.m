function [p2_umin, p2_umax, p2_M, p2_Np, p2_u_etfe, p2_y_etfe, p2_omega, p2_estimate_no_window, p2_estimate_windowed, p2_gamma_best] = HS2019_SysID_midterm_p2_sol_upload()
    % Sample solution of HS2019_SysID_midterm_p2.
    % Mingzhou Yin, 23 October 2019

    %% Given: plant function, validation input function
    sat_range = [-5,5];                         % search range for saturation
    sat_res = 0.1;                              % search resolution for saturation
    freq_res = pi/200;                          % required freq. resolution [rad/sample]
    var_w = 0.5^2;                              % variance of the Gaussian noise
    var_g = 0.2^2;                              % required variance of the estimate
    gamma_grid = 50:50:500;                     % search grid of gamma
    np_cv = 6;                                  % number of periods in validation data
    
    legi_num = 19951482;                        % your legi-number
    n_trans = 1;                                % number of periods to throw away in id data
    n_transcv = 1;                              % number of periods to throw away in cv data

    %% Question 1: input saturation
    n_test = 1e3;                                   % number of tests the find out the saturation
    sat_grid = sat_range(1):sat_res:sat_range(2);   % possible saturation values
    flag = zeros(size(sat_grid));                   % if the output is different up to some noise
    y0 = HS2019_SysID_midterm_p2_system_sim(legi_num,sat_grid(1)*ones(1,n_test));
    for i = 2:length(sat_grid)
        y = HS2019_SysID_midterm_p2_system_sim(legi_num,sat_grid(i)*ones(1,n_test));    % feed a step input to the system
        if abs(sum(y-y0))>4*sqrt(2*n_test*var_w)    % check if the output is significantly different to the previous one
            flag(i) = 1;
        end
        y0 = y;
    end
    idx = find(flag);
    u_min = sat_grid(idx(1))-sat_res;
    u_max = sat_grid(idx(end));

    %% Question 2: PRBS design
    M = 2^(ceil(log(2*pi/freq_res+1)/log(2)))-1;
    phi_u = ((u_max-u_min)/2)^2*(M+1)/M;            % input power spectral density of one period of PRBS (identical at all freq. except 0)
    np = ceil(var_w/phi_u/var_g);
    u_id = repmat(idinput(M,'prbs'),np+n_trans,1)*(u_max-u_min)/2+(u_max+u_min)/2;      % scale the PRBS input to [umin,umax]
    y_id = HS2019_SysID_midterm_p2_system_sim(legi_num,u_id);

    %% Question 3: unsmoothed estimate
    u_idt = u_id(n_trans*M+1:end);    	% truncate first n_trans periods
    y_idt = y_id(n_trans*M+1:end);
    
    % Method 1: frequency selection
    Y_N = fft(y_idt);
    U_N = fft(u_idt);
    idx_h = (0:M-1)*np+1;               % consider only harmonic frequencies w.r.t. the period
    Y_N = Y_N(idx_h);
    U_N = U_N(idx_h);
    G_N = Y_N./U_N;
    
%     % Method 2: time-domain averaging
%     u_avg = mean(reshape(u_idt,[],np),2);
%     y_avg = mean(reshape(y_idt,[],np),2);
%     Y_N = fft(y_avg);
%     U_N = fft(u_avg);
%     G_N = Y_N2./U_N2;
%     
%     % Method 3: frequency-domain averaging
%     u_trun = reshape(u_idt,[],np);
%     y_trun = reshape(y_idt,[],np);
%     Y_N = zeros(M,np);
%     U_N = zeros(M,np);
%     G_N = zeros(M,np);
%     for i = 1:np
%         Y_N(:,i) = fft(y_trun(:,i));
%         U_N(:,i) = fft(u_trun(:,i));
%         G_N(:,i) = Y_N(:,i)./U_N(:,i);
%     end
%     Y_N = mean(Y_N3,2);
%     U_N = mean(U_N3,2);
%     G_N = mean(G_N3,2);

    omega = 2*pi/M*(0:M-1)';
    idx = find(omega>0 & omega<pi);     % (0,pi)
    omega_p = omega(idx);
    G_N_p = G_N(idx);

    %% Question 4: smoothing + cross validation
    % a) smoothing
    Hann_fd = zeros(M,length(gamma_grid));
    G_N_smoothed = zeros(M,length(gamma_grid));
    U_N2 = abs(U_N).^2;                 % input power weighting
    for i = 1:length(gamma_grid)
        [~,Hann_fd(:,i)] = WfHann(gamma_grid(i),M);
        Hann_fd(:,i) = [Hann_fd(floor(M/2)+1:end,i);Hann_fd(1:floor(M/2),i)];
        
        % Calculate circular convolution manually
        for wn = 1:M
            Wnorm = 0;
            for xi = 1:M
                widx = mod(xi-wn,M)+1;
                G_N_smoothed(wn,i) = G_N_smoothed(wn,i) + Hann_fd(widx,i)*G_N(xi)*U_N2(xi);
                Wnorm = Wnorm + Hann_fd(widx,i)*U_N2(xi);
            end
            G_N_smoothed(wn,i) = G_N_smoothed(wn,i) / Wnorm;
        end
        
%         % Calculate circular convolution by cconv
%         G_N_smoothed(:,i) = (cconv(Hann_fd(:,i),abs(U_N).^2.*G_N,M))./(cconv(Hann_fd(:,i),abs(U_N).^2,M));
    end

    % b) cross validation
    u_cv = HS2019_SysID_midterm_p2_validation(M);               % load validation input
    y_cv = HS2019_SysID_midterm_p2_system_sim(legi_num,u_cv);
    u_cvt = u_cv(n_transcv*M+1:end);                         	% truncate out first n_transcv periods
    y_cvt = y_cv(n_transcv*M+1:end);
    Y_cv = fft(y_cvt);
    U_cv = fft(u_cvt);
    idx_hcv = (0:M-1)*(np_cv-n_transcv)+1;                      % consider only harmonic frequencies w.r.t. the period
    Y_cv = Y_cv(idx_hcv);
    U_cv = U_cv(idx_hcv);
    
    idx_cv = find(omega>0 & omega<pi & abs(U_cv)>max(abs(U_cv))/100);       % only do cross validation for freq. with enough input power

    error = zeros(length(gamma_grid),1);                        % error is defined as ||Y-\tilde{G}*U||_2 as in the problem
    for i = 1:length(gamma_grid)
        for j = idx_cv'
            error(i) = error(i)+abs(Y_cv(j)-G_N_smoothed(j,i)*U_cv(j))^2;
        end
        error(i) = sqrt(error(i));
    end
    [~,idx_gamma] = min(error);
    gamma_best = gamma_grid(idx_gamma);
    figure(11)
    plot(error)
    grid on; hold on;
    xlabel('\gamma')
    ylabel('Error');
    
	%% Comparison    
    figure(12)
    loglog(omega_p,abs(G_N_p));
    hold on
    loglog(omega_p,abs(G_N_smoothed(idx,idx_gamma)));
    
    xlim([omega_p(1),omega_p(end)]);
    ylim([1e-2,1e1])
    xlabel('Frequency (rad/sample)')
    ylabel('Magnitude');
    legend('Unsmoothed estimate','Optimal smoothed estimate')

    %% Outputs
    p2_M = M;
    p2_Np = np;
    p2_u_etfe = u_id;
    p2_umin = u_min;
    p2_umax = u_max;
    p2_y_etfe = y_id;
    p2_omega = omega_p;
    p2_estimate_no_window = G_N_p;
    p2_estimate_windowed = G_N_smoothed(idx,idx_gamma);
    p2_gamma_best = gamma_best;
end