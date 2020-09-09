a = 0.5;
b = 1;
% N = 4096;
p=2;


uniformnoise = false;

% u = randn(N, 1);    % i-th element represent u(i-1)
% if uniformnoise
%     w = (rand(N,1)-0.5) * 2 * sqrt(0.6);    % i-th element represent w(i)
% else
%     w = sqrt(0.2) * randn(N,1);
% end

% phi = zeros(N,1);
% phi(1,1) = 0;   % y(0)
% phi(:,2) = u;
% y = zeros(N,1);
no_n =1;
for N = [64 256 1024 4096 16384]


    no_try = 1000;
    theta_est = zeros(no_try, 2);
    rss_error = zeros(no_try, 1);
    for iter = 1:no_try

        u = randn(N, 1);    % i-th element represent u(i-1)

        if uniformnoise
            w = (rand(N,1)-0.5) * 2 * sqrt(0.6);    % i-th element represent w(i)
        else
            w = sqrt(0.2) * randn(N,1);
        end

        phi = zeros(N,2);
        phi(1,1) = 0;   % y(0)
        phi(:,2) = u;
        y = zeros(N,1);

        for i = 1:N
            y(i) = a*phi(i,1) + b*phi(i,2) + w(i);
            if i < N
                phi(i+1) = y(i);
            end
        end
        theta_est(iter,:) = phi\y;
        rss_error_total(iter) = sum((y-phi*theta_est(iter,:)').*(y-phi*theta_est(iter,:)'));
    end
    figure(no_n*2-1)
    subplot(2,1,1)
    hist(theta_est(:,1), linspace(0.3, 0.7, 50))
    subplot(2,1,2)
    hist(theta_est(:,2), linspace(0.8, 1.2, 50))

    figure(no_n*2)
    subplot(3,1,1)
    hist(rss_error_total/0.2);
    subplot(3,1,2)
    hist(chi2rnd(N-p, no_try,1));
    subplot(3,1,3)
    hist((N-p)*ones(no_try,1)+sqrt(2*(N-p))*randn(no_try,1));
    disp(N)



    no_n = no_n + 1;


end