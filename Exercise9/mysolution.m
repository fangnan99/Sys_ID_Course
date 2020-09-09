clear all
clc

load('Data_ex9.mat');

u = ex9_u;
y = ex9_y;
k = length(u);

[theta_e, fval] = fmincon(@(theta_e)PLRObjective(theta_e), zeros(k+4,1), [], [], [], [], [], [],...
                         @(theta_e)PLRConstraint(u, y, theta_e));
%% 1




%----------------------------------------------------------------------%
%% 2

% Starting from Linear Regression
k = length(u);

% One way to calculate theta_ls, using fmincon
% [theta_e_ls, fval] = fmincon(@(theta_e_ls)LRObjective(theta_e_ls), zeros(k+3,1), [], [], [], [], [], [],...
                        %  @(theta_e_ls)LRConstraint(u, y, theta_e_ls));
% theta_ls = theta_e_ls(1:3);


% Another way to calculate theta_ls                         
varphi = zeros(k, 3);
varphi(1,:) = [0 0 0];
varphi(2,:) = [-y(1) 0 u(1)];
for i=3:k
    varphi(i,:) = [-y(i-1) -y(i-2) u(i-1)];
end
theta_ls = varphi\y;

A_cap = [1 theta_ls(1), theta_ls(2)];
B_cap = [theta_ls(3), 0];

G_cap = tf(B_cap, A_cap, -1);

x = lsim(G_cap, u);

varzeta = zeros(k, 3);
varzeta(1,:) = [0 0 0];
varzeta(2,:) = [-x(1) 0 u(1)];
for i=3:k
    varzeta(i,:) = [-x(i-1) -x(i-2) u(i-1)];
end

theta_iv = (varzeta'*varphi)\(varzeta'*y);



















%----------------------------------------------------------------------%
%% 1.
% varphi(k) = [-w(k-1) -w(k-2) u(k-1) e(k-1)]
% where w(k) = -a1*w(k-1) - a2*w(k-2) + b1*u(k-1)
function f = PLRObjective(theta_e)  % theta = theta_e(1:4), e = theta_e(5:end)
% 2-norm of e
e = theta_e(5:end);
f = sqrt(e'*e);    
end

function [c, ceq] = PLRConstraint(u, y, theta_e)
c = [];
a_1 = theta_e(1);
a_2 = theta_e(2);
b_1 = theta_e(3);
c_1 = theta_e(4);
theta = theta_e(1:4);
e = theta_e(5:end);
k = length(u);
w = zeros(k, 1);
w(1) = 0;
w(2) = b_1 * u(1);
for i=3:k
    w(i) = -a_1*w(i-1) - a_2*w(i-2) + b_1*u(i-1);
end

varphi = zeros(k, 4);
varphi(1,:) = [0 0 0 0];
varphi(2,:) = [0 0 u(1) e(1)];
for i = 3:k
    varphi(i,:) = [-w(i-1) -w(i-2) u(i-1) e(i-1)];
end

% y_cap = zeros(k,1);
% y_cap(1) = w(1);
% for i=2:k
%     y_cap(i) = w(i) + c_1*e(i-1);
% end

y_cap = varphi*theta;
ceq = y - y_cap - e;
end









%----------------------------------------------------------------------%
%% 2.
function f = LRObjective(theta_e_ls)  % theta = theta_e_ls(1:3), e = theta_e_ls(4:end)
% 2-norm of e
e = theta_e_ls(4:end);
f = sqrt(e'*e);    
end

function [c, ceq] = LRConstraint(u, y, theta_e_ls)
c = [];

theta = theta_e_ls(1:3);
a_1 = theta(1);
a_2 = theta(2);
b_1 = theta(3);

e = theta_e_ls(4:end);
k = length(u);
% w = y_cap
w = zeros(k, 1);
w(1) = 0;
w(2) = -a_1*y(1) + b_1 * u(1);
for i=3:k
    w(i) = -a_1*y(i-1) - a_2*y(i-2) + b_1*u(i-1);
end

ceq = y - w - e;
end