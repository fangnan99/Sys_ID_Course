function [p3_b_ML,p3_b_MAP,p3_cv_error,p3_prior_best] = HS2019_SysID_final_p3_11111111()

%% Solution for Problem 3

%% General instructions for solution
% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the null variables as already
% provided in the code

% Extract Legi from Filename
name = mfilename;
LegiNumber = str2double(name(end-7:end));

% Obtain experiment data
[p3_u,p3_y,p3_u_cv,p3_y_cv] = HS2019_SysID_final_p3_GenerateData(LegiNumber);

%% Task 1: Maximum likelihood estimate





p3_b_ML         = zeros(8,1);   % vector of dimension 8x1


%% Task 2: Maximum a posteriori estimates





p3_b_MAP        = zeros(8,5);   % matrix of dimension 8x5
p3_cv_error    	= zeros(5,1);   % vector of dimension 5x1
p3_prior_best   = 0;            % scalar integer in the set {1,2,3,4,5}


disp('************************************************************')
disp('**                       Problem 3                        **')
disp('**                 --------------------                   **')
disp('************************************************************')

%% Define Matrices
N = length(p3_u);
dim_b = 8;
sigma_e = 0.5^2 * eye(N);
inv_sigma_e = inv(sigma_e);
Mu_e = zeros(N, 1);

phi = zeros(N, 8);
for i = 1:N
    if i < 9
        for j = 1:i-1
            phi(i, j) = p3_u(i-j);
        end
    else
        for j = 1:8
            phi(i,j) = p3_u(i-j);
        end
    end
end

%% Part 1
mytext=sprintf('\n\n\nPart a\n');
disp(mytext)

b_est_ml = (phi'*inv_sigma_e*phi) \ (phi'*inv_sigma_e*p3_y);

mytext=sprintf(['Let phi(k) = [u(k-1) u(k-2) ...u(k-8)], then y(k) = phi(k)*b + e(k).\n',...
                'Let Y=[y(0); y(1); ...; y(99)], Phi=[phi(0); phi(1); ...; phi(99)], E=[e(0); e(1); ...; e(99)],\n',...
                'then Y = Phi*b + E.\n\n',...
                'The maximun likelihood estimate of b is defined as\n',...
                'b_hat_ML = argmax(b){P(Y|b)}\n',...
                'e(k) is i.i.d. and e(k) ~ N(0, sigma^2), so E ~ N(0, Sigma=sigma^2*I).\n',...
                'Therefore, Y ~ N(Phi*b, Sigma), the likelihood function is\n',...
                'P(Y|b) = 1/sqrt((2*pi)^100*det(Sigma)) * exp(-1/2 * (Y-Phi*b)'' * Sigma^(-1) * (Y-Phi*b))\n',...
                'By simple algebra, we have\n',...
                'b_hat_ML = argmax(b){P(Y|b)} = argmin(b){b''*Phi''*Sigma^(-1)*Phi*b - 2Y''Sigma^(-1)*Phi*b}\n',...
                'According to the hints, b_hat_ML is the solution to the following equation about b\n',...
                'Phi'' * Sigma^(-1) * Phi * b - Phi'' * Sigma^(-1) * Y = 0\n',...
                'Thus, b_hat_ML = (Phi'' * Sigma^(-1) * Phi) \\ (Phi'' * Sigma^(-1) * Y).\n',...
                ]);
disp(mytext)

%% Part 2
mytext=sprintf('\n\n\nPart b\n');
disp(mytext)

b_est_map = zeros(dim_b, 5);
S = zeros(5, dim_b, dim_b);
S(1,:,:) = eye(dim_b);
for i = 1:dim_b
    S(2,i,i) = 0.8^i;
    S(3,i,i) = 0.5^i;
end
for i = 1:dim_b
    for j = 1:dim_b
        S(4,i,j) = 0.8^(max(i,j));
        S(5,i,j) = 0.5^(max(i,j));
    end
end

for i = 1:5
    b_est_map(:,i) = (inv(reshape(S(i,:,:),dim_b,dim_b)) + phi'*inv_sigma_e*phi) \ (phi'*inv_sigma_e*p3_y);
end

mytext=sprintf(['The maximun a posteriori estimate of b is defined as\n',...
                'b_hat_MAP = argmax(b){P(b|Y)}\n',...
                'According to Bayes'' theorem,\n',...
                'b_hat_MAP = argmax(b){P(b|Y)} = argmax(b){P(Y|b)P(b)/P(Y)} = argmax(b){P(Y|b)P(b)}',...
                'and\n',...
                'P(Y|b)P(b) = 1/sqrt((2*pi)^100*det(Sigma))*exp(-1/2*(Y-Phi*b)''*Sigma^(-1)*(Y-Phi*b))*1/sqrt((2*pi)^8*det(S))*exp(-1/2*b''*S^(-1)*b)\n',...
                'By simple algebra, the optimization problem could be written as\n',...
                'b_hat_MAP = argmax(b){P(Y|b)P(b)} = argmin(b){b''*S^(-1)*b + b''*Phi''*Sigma^(-1)*Phi*b - 2Y''Sigma^(-1)*Phi*b}\n',...
                'According to the hints, b_hat_MAP is the solmution to\n',...
                '[S^(-1) + Phi'' * Sigma^(-1) * Phi] * b - Phi'' * Sigma^(-1) * Y = 0\n',...
                'Thus, b_hat_MAP = (S^(-1) + Phi'' * Sigma^(-1) * Phi) \\ (Phi'' * Sigma^(-1) * Y).\n',...
                ]);
disp(mytext)

%% Part 3
epsilon = zeros(5,1);
N_cv = length(p3_u_cv);
phi_cv = zeros(N, 8);
for i = 1:N_cv
    if i < 9
        for j = 1:i-1
            phi_cv(i, j) = p3_u_cv(i-j);
        end
    else
        for j = 1:8
            phi_cv(i,j) = p3_u_cv(i-j);
        end
    end
end



% figure; hold on;
% plot(p3_y_cv, 'color', 'k', 'LineWidth', 2)

% Debug Code
%% See the cv error of ML estimator
% y_pred_ml = phi_cv * b_est_ml;
% y_pred_ml_error = p3_y_cv - y_pred_ml;
% plot(y_pred_ml)
% disp(mean(y_pred_ml_error.*y_pred_ml_error))

for i = 1:5
    y_pred = phi_cv * b_est_map(:,i);
    % plot(y_pred)
    y_pred_error = p3_y_cv - y_pred;
    epsilon(i) = mean(y_pred_error.*y_pred_error);
end
% legend('true', 'ml', '1', '2', '3', '4', '5');
best_idx = find(epsilon == min(epsilon));

mytext=sprintf(['\nDifference between ML and MAP\n',...
                'The difference of ML and MAP estimate of b is that the optimization objective functions\n',...
                'are different. In ML estimate the likelihood function P(Data|Parameter) is maximized, but\n',...
                'in MAP estimate the posteriori function P(Parameter|Data) is maximized, or (because of \n',...
                'Bayes''s rule)equivalently the product of priori and likelihood is maximized. Therefore\n',...
                'the ML estimate is purely based on data observed while the MAP estimate also uses prior\n',...
                'knowledge or assumptions about the parameters to regularize the optimization. In case that\n',...
                'the prior distribution is uniform the two estimates would be exactly the same. In this\n',...
                'particular problem, b and E are both normally distributed, the optimization objective\n',...
                'function of MAP estimate has and an additional term, b''*S^(-1)*b, which causes the two\n',...
                'estimates to be different.\n',...
                ]);
disp(mytext)

% plot 2 Gaussion distributions to show effect of S
mu_b = [0;0];
cov1_b = [S(1,1,1) S(1,1,7); S(1,1,7) S(1,7,7)];
cov3_b = [S(3,1,1) S(3,1,7); S(3,1,7) S(3,7,7)];

points_x = [b_est_ml(1); b_est_map(1,1); b_est_map(1,3)];
points_y = [b_est_ml(7); b_est_map(7,1); b_est_map(7,3)];

inv_sigma3 = inv(cov3_b);
inv_sigma1 = inv(cov1_b);
[X,Y]=meshgrid(-1.5:0.1:1.5,-1.5:0.1:1.5);

a = exp(-1/2 * (inv_sigma3(1,1)*X.*X + inv_sigma3(1,2)*X.*Y + inv_sigma3(2,1)*Y.*X + inv_sigma3(2,2)*Y.*Y));
P_prior_3= 1/sqrt((2*pi)^2 * det(cov3_b)) * a;
P_prior_3=reshape(P_prior_3,size(X));
a = exp(-1/2 * (inv_sigma1(1,1)*X.*X + inv_sigma1(1,2)*X.*Y + inv_sigma1(2,1)*Y.*X + inv_sigma1(2,2)*Y.*Y));
P_prior_1= 1/sqrt((2*pi)^2 * det(cov1_b)) * a;
P_prior_1=reshape(P_prior_1,size(X));

N_plot_points = length(points_x);
points_p1 = zeros(N_plot_points,1);
points_p3 = zeros(N_plot_points,1);
for i=1:N_plot_points
    points_p1(i) = 1/sqrt((2*pi)^2 * det(cov1_b)) * exp(-1/2 * ([points_x(i); points_y(i)]-mu_b)'*inv_sigma1*([points_x(i); points_y(i)]-mu_b));
    points_p3(i) = 1/sqrt((2*pi)^2 * det(cov3_b)) * exp(-1/2 * ([points_x(i); points_y(i)]-mu_b)'*inv_sigma3*([points_x(i); points_y(i)]-mu_b));
end
figure(1);
hold on;
mesh(X,Y,P_prior_1)
plot3(points_x(1), points_y(1), points_p1(1), '*', 'MarkerSize', 8, 'color', 'b')
plot3(points_x(2), points_y(2), points_p1(2), 'o', 'MarkerSize', 8, 'color', 'r')
xlabel('b_1');
ylabel('b_7');
zlabel('Prior Probability');
title('Prior Probability based on S_1');
legend('Prior PDF','ML estimate', 'MAP estimate')
view(70, 45);

figure(2);
hold on;
mesh(X,Y,P_prior_3)
plot3(points_x(1), points_y(1), points_p3(1), '*', 'MarkerSize', 8, 'color', 'b')
plot3(points_x(3), points_y(3), points_p3(2), 'o', 'MarkerSize', 8, 'color', 'r')
xlabel('b_1');
ylabel('b_7');
zlabel('Prior Probability');
title('Prior Probability based on S_3');
legend('Prior PDF','ML estimate', 'MAP estimate')
view(70, 45);

mytext=sprintf(['\nEffect of choices of S\n',...
                'In this problem, S is the prior assumption of b''s covarince, different choice of S would affect\n',...
                'the prior distribution of b, thus affect the MAP estimate of b. Because ML estimate maximizes\n',...
                'P(Y|b) and MAP estimate maximizes P(Y|b)P(b), the MAP estimate of vector b tend to make P(b) larger\n',...
                'comparing to the ML estimate, which means the MAP estimate of elements in b tend to have more\n',...
                'similar correlations as we assumed in S.\n',...
                'Figure 1 and Figure 2 the prior distribution of b_1 and b_7 based on S1 and S3 respectively. In the\n',...
                'plots, the ML and corresponding MAP estmates are denoted with markers. We can see thatwhen S=S1 the\n',...
                'estimates are close to each other. Since the prior distribution is relatively "flat", it does not\n',...
                'have much influence on the posteriori. However, when S=S3 the MAP estimate is further away\n',...
                'from the ML estimator because the prior distribution is "sharp", which means we have strong prior\n',...
                'assumptions about the distribution of b.\n',...
                'In the validation, the MAP estimate based on S_2 has the smallest MSE, which means it is a\n',...
                'relatively better assumption of the covarince of b. On the contrary, S_3 gives the largest\n',...
                'predicion MSE on the cross validation data, which suggests it might be a bad choice and may bring\n',...
                'bias to the estimate.\n']);
disp(mytext)

%% Return values
p3_b_ML = b_est_ml;
p3_b_MAP = b_est_map;
p3_cv_error = epsilon;
p3_prior_best = best_idx;

end



