% This function is made by Ahmed Tahan
%{
    
    It's intended to apply the self-tuning regulator for a given system
    such as 
                                            y         z^(-d) Bsys
                                Gp = ------ = ----------------------
                                            u               Asys

    the controller is given in the form of 
                                    
                                            u          S 
                                Gc = ------ = ------- 
                                           err         R

    the closed loop transfer function
           y                 z^(-d)BsysS                    z^(-d)BsysS          z^(-d)BsysS
        ------ = ---------------------------------- =  -------------------  = ------------------------
          uc             AsysR + z^(-d)BsysS           Am A0                   alpha

where 
-- y : output of the system
-- u : control action (input to the system)
-- uc : required output (closed loop input-reference, command signal)
-- err = error between the required and the output --> = uc - y
-- Asys = 1 + a_1 z^-1 + a_2 z^-1 + ... + a_na z^(-na) 
-- Bsys = b_0 + b_1 z^-1 + b_2 z^-1 + ... + b_nb z^(-nb)
-- R = 1 + r_1 z^-1 + r_2 z^-1 + ... + r_nr z^(-nr) --> [1,  r_1,  r_2,  r_3,  ..., r_nr] 
-- S = s_0 + s_1 z^-1 + s_2 z^-1 + ... + s_ns z^(-ns) --> [s_0,  s_1, s _2,  s_3,  ..., s_ns] 
-- d : delay in the system. Notice that this form of the Diaphontaing solution
        is available for systems with d>=1
-- Am = required polynomial of the model = 1+m_1 z^-1 + m_2 z^-1 + ... + m_nm z^(-m_nm)
-- A0 = observer polynomail for compensation of the order = 1 + o_1 z^-1 + o_2 z^-1 + ... + o_no z^(-no)
-- alpha:required characteristic polynomial = Am A0 = 1 + alpha1 z^-1 + alpha2 z^-1 + ... + alpha_(nalpha z)^(-nalpha) 

Steps of solution:
1- initialization of the some parameters (theta0, P, Asys, Bsys, S, R, T, y, u, err, dc_gain).
2- assume at first the controller is unity. Get u, y of the system
3- RLS and get A, B estimated for the system. 
4- Solve the Diophantine equation using A, B and the specified "alpha = AmA0" and get S, R of the controller.
5- find "u" due to this new controller and then y. 
6- repeat from 3 till the system converges.

Function Inputs and Outputs
Inputs
    uc: command signal (column vector)
    Asys = [1,  a_1,  a_2,  a_3,  ..., a_na] ----> size(1, na)
    Bsys = [b_0,  b_1, b _2,  b_3,  ..., a_nb]----> size(1, nb)
    d : delay in the system (d>=1)
    Ts : sample time (sec.)
    Am = [1,  m_1,  m_2,  m_3,  ..., m_nm]---> size(1, nm) 
    A0 = [1,  o_1,  o_2,  o_3,  ..., o_no]---> size(1, no) 
    alpha : [1,  alpha_1, alpha _2,  alpha_3,  ..., alpha_nalpha] ----> size(1, nalpha) row vector

Outputs
Theta_final : final estimated parameters
Gz_estm : estimated pulse transfer function
Gc = the controller by Diophantine equation = S/R
Gcl = closed loop transfer function 

here we may input the signal uc by dividing it by the dc gain in order to
force the output to follow the input in magnitude

%}



function [ Theta_final Gz_estm, Gc, Gcl] = ISTR( uc, Asys, Bsys, d, Ts, Am, A0 )

na = length(Asys) - 1;
nb = length(Bsys) - 1;
nu = na + nb + 1; % number of unkowns
nr = nb+d-1; 
ns = na-1;
nalpha = na + nb +d -1; % systems final order
N = length(uc); % length of the vectors
n = max(na+1, nb+d+1); % used to avoid zeros at the begining and added 1 to avoid zero index in matlab
final_time = (N-1)*Ts;
t = 0:Ts:final_time;

%% Initializations
% Covariance matrix Initialization
alph = 10E6;
for i = 1:(n-1)
P{1,i} = alph*eye(nu, nu);
end

% Initialization of estimated parameters
theta = cell(1, N);
theta_0 = 0.1;
A = 1;
for i = 1:(n-1)
%     theta{1, i} = zeros(nu,1);
    theta{1, i} = theta_0*ones(nu,1);   % initial estimation should be a small number but not zero 
        %in order to have at first S and R that are reasonable and you have to 
        % check if that they will converge to the real paramaters or not
    for k = 1: na
        A(k+1) = theta{1, i}(k);
    end
        
    for k = na+1:nu
        B(k-na) = theta{1, i}(k);
    end
end

% Initialization of the system
u = zeros(length(uc), 1);
y_estm = zeros(length(uc), 1);
err = zeros(length(uc), 1);
y_final = zeros(length(uc), 1);

% Initialization of R, S
R = 1;
S = 1;
for i = 1:nr
    R(i+1) = 0;
end
    
for i = 1:ns
    S(i+1) = 0;
end

% initial dc gain
dc_gain = sum(B)/sum(A);

%closed loop initilization by taking into consideration the term z^(-d)*B*S in the denominator by the term DBS_sys
X = conv(Bsys, S);
DBS_sys = zeros(length(X)+d, 1);
DBS_sys(d+1:length(X)+d) = X';
Acl = conv(Asys, R)+ DBS_sys';
Bcl = DBS_sys';

% Initialization of cell for legend the estimated paramters
Name = cell(nu, 1);

% Required Characterisitic Polynomial
alpha = conv(Am, A0);

%% Recursive Least Squares and Diophantine Controller Estimation Algorithm
for j = n : N
    X = conv(Bsys, S);
    DBS_sys = zeros(length(X)+d, 1);
    DBS_sys(d+1:length(X)+d) = X';
    Acl = conv(Asys, R)+ DBS_sys';
    Bcl = DBS_sys';
    y_final(j) = outputestimation( Acl, Bcl/dc_gain, d, uc, y_final, j );
        
    % Calculation of the dc gain of the closed loop tf based on the estimated A, B, S, R
    X = conv(B, S);
    DBS = zeros(length(X)+d, 1);
    DBS(d+1:length(X)+d) = X';
    num_gain = sum(conv(B, S));
    den_gain = sum(conv(A, R)+ DBS');
    dc_gain = num_gain/den_gain;
               

    % Closed Loop Estimation Iteration [ loop estimation (assume at first u = 0, for the error --> err(j) = uc(j) - y(j-1), y(j-1)
    % is calculated based on the previous input and then find the new action u(j) ]   
    y_estm(j) = outputestimation( Asys, Bsys, d, u, y_estm, j );
    err(j) = uc(j)/dc_gain - y_estm(j); % dc gain must be added to the uc "commanded signal" as the dc gain is defined between the clolsed transfer functoin between y and uc
    u(j)  = outputestimation( R, S, 0, err, u, j );                
    k(j) = y_estm(j)/dc_gain;

    % Filling the Phi matrix
    for i = 1 : na % this for loop used to fill parts in the same line that pertains to the output "denomenator"
        if ((j-i)<=0)
            Phi(j, i) = 0;
        else
            Phi(j, i) = -y_estm((j-i));
        end
    end
    
        
    for i = na+1 : nu % this for loop used to fill parts in the same line that pertains to the input "numerator" starts from the last postion column of "na +1 
        if ((j-(i-na)-d)<=0) % add na as we left the output going to the input and we start index after "na"
            Phi(j, i) = 0;
        else
            Phi(j, i) = u((j-d-i + (na+1)));
        end
    end
        
    % Phi = [phi(i) phi(i+1) ... phi(N)], Phi is the matrix N*nu but the phi is the vector nu*1
    phi = Phi(j, :)'; % to get the j th row and then convert into column and that what we need in RLS not Phi
        
    % RLS algorithm
    P{1, j} = P{1, j-1} - P{1, j-1}*phi*inv(1+phi'*P{1, j-1}*phi)*phi'*P{1, j-1}; % Covariance Matrix
    K = P{1, j}*phi; % Gain
    theta{1, j} = theta{1, j-1} + K*(y_estm(j)-phi'*theta{1, j-1}); % estimated parameters
           

    %  1DOF Controller
     
    % Preparation of A, B for Diophantine Equation
    A = 1;
    for i = 1: na
        A(i+1) = theta{1, j}(i);
    end
        
    for i = na+1:nu
        B(i-na) = theta{1, j}(i);
    end
           
    [ S, R ] = Diophantine( A, B, d, alpha );
            
%             A
%             B
%             R
%             S

end
    
Theta = cell2mat(theta); % converting the Total_loop_estimate from cell array to matrix, the final column is the same as the final estimate
Theta_final = Theta(:, N); % final estimation of parameters

%% Parameters and transfer function preparation 

%Closed loop
X = conv(B, S);
DBS = zeros(length(X)+d, 1);
DBS(d+1:length(X)+d) = X';
Acl = conv(Asys, R)+ DBS';
Bcl = DBS'/dc_gain;
Gcl = tf(Bcl, Acl, Ts);

%System
z = tf('z');
 AA = z^(d+nb);
 BB = 0;
for i = 1:na
    a(i) = Theta_final(i);
    AA = AA + a(i)*z^(d-na+nb+na-i);
end
for i = na+1:nu
    b(i-na) = Theta_final(i);
    BB = BB + b(i-na)*z^(nb - (i-(na+1)));
end

Gz_estm = BB/AA;
[num, den] = tfdata(Gz_estm);
num = cell2mat(num);
den = cell2mat(den);
Gz_estm = tf(num, den, Ts);

%Controller
 z = tf('z');
 RR = z^(nr);
 SS = 0;
for i = 1:nr
    r(i) = R(i+1);
    RR = RR + r(i)*z^(nr-i);
end
for i = 1:(ns+1)
    s(i) = S(i);
    SS = SS +s(i)*z^(nr+1-i);
end
Gc = SS/RR;
[num, den] = tfdata(Gc);
num = cell2mat(num);
den = cell2mat(den);
Gc = tf(num, den, Ts);

%% Validation

% Output Estimated System
figure(2)
hold all
plot(t, y_estm, 'LineWidth', 2)
% plot(t, y_final, 'LineWidth', 2);
grid on 

% Required Output (Input uc)
plot(t, uc, 'LineWidth', 2);
legend(['Output Estimated System'], ['Required Output (Input uc) ',num2str(na), ' Order']);
xlabel('Time (sec.)');
ylabel('Output');
title('Output Vs Required Output')

% Parameter History Plot
for i = 1 : nu
    figure(3)
    hold all
    plot(t, Theta(i, :), 'LineWidth', 2);
    if (i>=1 & i<=na)
        Name{i, 1} = ['a', num2str(i)]; % Polynomal A parameters
    else
        Name{i, 1} = ['b', num2str(i-na-1)]; % Polynomal B parameters
    end
end

grid on
xlabel('Time (sec.)')
ylabel('Estimated Parameters')
title('Estimated Parameters History')
legend(Name)

figure(4)
plot(t, u, 'LineWidth', 2)
grid on
xlabel('Time (sec.)')
ylabel('Control action')
title('Control Action signal Vs Time')

print(2,'-dpng','-r500','Output')
print(3,'-dpng','-r500','Parameters')
print(4,'-dpng','-r500','Control_action')
end