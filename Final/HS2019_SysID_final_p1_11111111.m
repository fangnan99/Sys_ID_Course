function [p1_theta_est1,p1_Phi,p1_theta_est2,p1_y_pred] = HS2019_SysID_final_p1_11111111()
%% Solution for Problem 1
%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= str2num(name(end-7:end));

[p1_u, p1_y, p1_theta_hat, p1_u_past, p1_y_past,p1_pred_err] = HS2019_SysID_final_p1_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variable p1_U to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Task 1: Obtain initial estimate

p1_theta_est1 = [];
p1_Phi = [];

%% Task 2: Improve your estimate

p1_theta_est2 = [];


%% Task 3: Compute prediction

p1_y_pred = [];

disp('************************************************************')
disp('**                       Problem 1                        **')
disp('**                 --------------------                   **')
disp('************************************************************')

%% Compute Re
N = length(p1_y);
Re = zeros(N, N);
for i = 1:N
    for j = 1:N
        Re(i,j) = CovarinceEIJ(i,j);
    end
end
%% Compute Phi
phi = zeros(N, 6);
phi(1,:) = [0 0 0 0 0 0];
phi(2,:) = [-p1_y(1) 0 0 p1_u(1) 0 0];
phi(3,:) = [-p1_y(2) -p1_y(1) 0 p1_u(2) p1_u(1) 0];
for i = 4:N
    phi(i,:) = [-p1_y(i-1) -p1_y(i-2) -p1_y(i-3) p1_u(i-1) p1_u(i-2) p1_u(i-3)];
end

%% Compute Rv
Rv = zeros(N, N);
for i = 1:N
    for j = 1:N
        Rv(i,j) = CovarinceVIJ(i,j);
    end
end

mytext=sprintf('\n\n\nPart 1\n');
disp(mytext)
%% Estimate theta

Z_est1 = inv(Re)*phi*inv(phi'*inv(Re)*phi);
theta_est1_blue = Z_est1'*p1_y;
cov_theta_est1_blue = Z_est1'*Rv*Z_est1;
mytext=sprintf(['The system model could be written as A(z)y(k) = B(z)u(k)+C(z)e(k), where C(z)=1+c_1*z^-1+c_2z^-2.\n',...
                'Under the assumption that c_1 = c_2 = 0, C(z) = 1 and the model becomes equation error model.\n',...
                'From the covarince between e at different time steps, we know that e(k) is correlated noise. Thus\n',...
                'to obtain a linear estimate theta_hat_est1 of the parameter, the BLUE estimator can be used.\n',...
                'Then a row of the regressor phi(k) would be [-y(k-1), -y(k-2), -y(k-3), u(k-1), u(k-2), u(k-3)]\n',...
                'and the regressor Phi would be [phi(1); phi(2); ...; phi(N)]\n',...
                'Then a BLUE estimator could be formed in the following form.\n',...
                'theta_hat_est1 = Z''*Y\n',...
                'where Z = inv(R_e)*Phi*inv(Phi''*inv(R_e)*Phi), Y is the output vector, and R_e is the covarince\n',...
                'matrix of noise vector E = [e(1); e(2); ...; e(N)].\n\n',...
                'The standard statistical features of a BLUE estimator is the following form according to Slide 9.32\n',...
                'E{theta_hat_est1} = theta_0    (i.e. the estimator is unbiased)\n',...
                'cov{theta_hat_est1} = inv(Phi''*inv(R_e)*Phi)\n',...
                'However, in this problem, the statistical features of this estimator is slightly different. Because\n',...
                'the true noise contained in the data v(k)=C(z)e(k) is unknown, and the estimator is formed based on\n',...
                'the assumption that c_1=c_2=0. Similar to the derivation of statistical properties for the standard\n',...
                'BLUE estimator, under the assumption that Phi is fixed and V is a random vector, we could derive that\n',...
                'this is still an unbiased estimator because V has zero mean, i.e.\n',...
                'E{theta_hat_est1-theta_0} = E{Z''(Phi*theta_0 + V) - theta_0} = E{theta_0 + Z''*V - theta_0} = 0.\n',...
                'However, the covarince of theta_hat_est1 is different, we could derive that\n',...
                'cov{theta_hat_est1} = Z''*R_v*Z    where R_v is the covarince matrix of filtered noise v(k)=C(z)e(k).\n']);
disp(mytext)

mytext=sprintf('\n\n\nPart 2\n');
disp(mytext)
%% Re-estimate theta

Z_est2 = inv(Rv)*phi*inv(phi'*inv(Rv)*phi);
theta_est2_blue = Z_est2'*p1_y;
cov_theta_est2_blue = Z_est2'*Rv*Z_est2;
mytext=sprintf(['Now that we know the values of c_1 and c_2, the estimator in Part 1 could be improved in either\n',...
                'of the following approaches:\n',...
                '1. Note that C(z) in stably invertible under given values of c_1 and c_2. Then we could filter\n',...
                '   y and u with 1/C(z) and get y_f and u_f. Then the model becomes y_f(k) = B(z)/A(z)u_f(k) + \n',...
                '   1/A(z)e(k). Then use the same approach as the previous problem and filtered data to get a new\n',...
                '   estimate.\n',...
                '2. Make use of c_1 and c_2 by calculating the covarince of filtered noise v(k) = C(z)e(k), and\n',...
                '   then use the same approach in the previous problem but replacing R_e with R_v.\n',...
                'It was varified that both approaches gives the same result.\n',...
                'Take the second method, then the new estimate is\n',...
                'theta_hat_est2 = Z''*Y\n',...
                'where Z = inv(R_v)*Phi*inv(Phi''*inv(R_v)*Phi), where R_v is covarince of filtered noise vector V.\n',...
                'R_v can be calculated because we know the accurate form of C(z) and covarince of noise vector E.\n\n',...
                'As for the statistical properties of the new estimator, if follows the standard form of BLUE estimator\n',...
                'because this estimator is a "true" BLUE estimator with no assumptions of the system. Therefore\n',...
                'E{theta_hat_est2} = theta_0    (i.e. the estimator is unbiased)\n',...
                'cov{theta_hat_est2} = inv(Phi''*inv(R_v)*Phi)\n',...
                'The mean of the new estimator is still unbiased, and its variance is not bigger than that of the\n',...
                'first one. The reason is that the new estimor is the "true" BLUE estimator and the previous one\n',...
                'is essentially a different linear unbiased estimator with a similar form. Thus by the property\n',...
                'of BLUE estimator, the new one has lower or equal covarince.\n']);
disp(mytext)

%% Another method: filter y and u with C^-1     ****** Result the same as calculating R_v ******
% C = [1 1 pi/4];
% G_inv_C = tf([1 0 0], C, -1);
% y_filtered = lsim(G_inv_C, p1_y);
% u_filtered = lsim(G_inv_C, p1_u);
% phi_filtered = zeros(N, 6);
% phi_filtered(1,:) = [0 0 0 0 0 0];
% phi_filtered(2,:) = [-y_filtered(1) 0 0 u_filtered(1) 0 0];
% phi_filtered(3,:) = [-y_filtered(2) -y_filtered(1) 0 u_filtered(2) u_filtered(1) 0];
% for i = 4:N
%     phi_filtered(i,:) = [-y_filtered(i-1) -y_filtered(i-2) -y_filtered(i-3) u_filtered(i-1) u_filtered(i-2) u_filtered(i-3)];
% end
% Z_est2_filtered = inv(Re)*phi_filtered*inv(phi_filtered'*inv(Re)*phi_filtered);
% theta_est_blue_filtered = Z_est2_filtered'*y_filtered;


%% Part 3 One Step Ahead Prediction
mytext = sprintf('\n\n\nPart 3\n');
disp(mytext)
% Bad code, wrong
% phi_prdt = [-p1_y_past(1) -p1_y_past(2) -p1_y_past(3) p1_u_past(1) p1_u_past(2) p1_u_past(3)];
% y_prdt = phi_prdt * p1_theta_hat;

% Based on Slide 10.18
osa_predictor = [p1_theta_hat; 1; pi/4];    % Coefficients [a1 a2 a3 b1 b2 b3 c1 c2]
phi_prdt = [-p1_y_past(1) -p1_y_past(2) -p1_y_past(3) p1_u_past(1) p1_u_past(2) p1_u_past(3) p1_pred_err(1) p1_pred_err(2)];
y_prdt = phi_prdt * osa_predictor;
mytext=sprintf(['In this part, a one-step-ahead prediction of ARMAX model is performed based on the equation on\n',...
               'Slide 10.18.\n',...
               'y_hat(k|theta) = B(z)u(k) + (1-A(z))y(k) + (C(z)-1)epsilon(k)\n',...
               'where epsilon(k) is the prediction error at time step k.\n',...
               'In this case, the u and y at past 3 time steps and prediction errors of past 2 time steps were\n',...
               'used in the calculation, and\n',...
               'y_pred(40) = [a1 a2 a3 b1 b2 b3 c1 c2] * [-y(39) -y(38) -y(37) u(39) u(38) u(37) epsilon(39) epsilon(38)]''\n']);
disp(mytext)

%% Return values
p1_theta_est1 = theta_est1_blue;
p1_Phi = phi;
p1_theta_est2 = theta_est2_blue;
p1_y_pred = y_prdt;
end


%% Functions
function covij = CovarinceEIJ(i, j)
    if i == j
        covij = 0.9;
    elseif abs(i-j) <= 3
        covij = 0.2;
    else
        covij = 0;
    end
end

function covij = CovarinceVIJ(i, j)
    vi_components = [i i-1 i-2];
    vi_coefficients = [1 1 pi/4];
    idx = find(vi_components > 0);
    vi_components = vi_components(idx);
    vi_coefficients = vi_coefficients(idx);
    
    vj_components = [j j-1 j-2];
    vj_coefficients = [1 1 pi/4];
    idx = find(vj_components > 0);
    vj_components = vj_components(idx);
    vj_coefficients = vj_coefficients(idx);
    
    covij = 0;
    for m = 1:length(vi_components)
        for n = 1:length(vj_components)
            covij = covij + CovarinceEIJ(vi_components(m), vj_components(n)) * vi_coefficients(m) * vj_coefficients(n);
        end
    end
    
end
