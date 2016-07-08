%   This function is made by Ahmed ElTahan. 
%{
It's intended to solve the Diphantine equation in the form of 

                            AR + z^(d) BS = A0Am = alpha;
where

-- A = 1 + a_1 z^-1 + a_2 z^-1 + ... + a_na z^(-na) --> [1,  a_1,  a_2,  a_3,  ..., a_na] 
-- B = b_0 + b_1 z^-1 + b_2 z^-1 + ... + b_nb z^(-nb) --> [b_0,  b_1, b _2,  b_3,  ..., a_nb] 
-- R = 1 + r_1 z^-1 + r_2 z^-1 + ... + r_nr z^(-nr) --> [1,  r_1,  r_2,  r_3,  ..., r_nr] 
-- S = s_0 + s_1 z^-1 + s_2 z^-1 + ... + s_ns z^(-ns) --> [s_0,  s_1, s _2,  s_3,  ..., s_ns] 
-- d : delay in the system. Notice that this form of the Diaphontaing solution
        is available for systems with d>=1

-- alpha = 1 + alpha1 z^-1 + alpha2 z^-1 + ... + alpha_(nalpha z)^(-nalpha)
        --> [1,  alpha_1, alpha _2,  alpha_3,  ..., alpha_nalpha]  required characteristic polynomial
-- Am = required polynomial of the model;
-- A0 = observer polynomail for compensation of the order

The function input outputs are given in the following

                        function [ S, R ] = Diophantine( A, B, d, alpha )

the functions is used to estimate the polynomials S and R which are the
numerator and the denomenator of the controller transfer function,
respectively.

The Solution is given in matrix form by solving a linear system of
equations such as 

                    M*theta = (V-Y) --> theta = M^(-1)*(V-Y)

-- M : Sylvester matrix
-- V: vector contains the "alpha" polynomail coefficients without "1" at the
        first of it. 
        
                            V = transpose([alpha_1, alpha _2,  alpha_3,  ..., alpha_nalpha])
                            size(nalpha, 1)

-- Y: vector contains the "A" polynomail coefficients without "1" at the
        first of it.

                            Y = transpose([a_1,  a_2,  a_3,  ..., a_na, 0, 0, ..., 0])
                            size(nalpha, 1)

-- theta: vector contains the unknowns. That is, the coefficients of the R
        polynomial and the coefficients of the S polynomial
    
                            theta = tranpose([r_1, r_2, ..., r_nr, s_0, s1, s2, ..., s_ns])

 %}



function [ S, R ] = Diophantine( A, B, d, alpha )
%% Calculating orders

if(d<=0)
    disp('error, the delay d must be >=1')
end
na = length(A)-1;
nb = length(B)-1;
nr = nb + d - 1;
ns = na - 1;
n_alpha = na + nb +d - 1;

%% initialization of M matrix and V matrix
M = zeros(n_alpha, n_alpha);
Y = zeros(n_alpha, 1);
V = zeros(n_alpha, 1);
Y(1:na, 1) = A(2:end);
V(1:n_alpha, 1) = alpha(2:end);

%% M matrix filling
for i = 1: n_alpha
    if(i<=nr)
        M(i:(length(A)+i-1), i) = A;
    else
        M(d+(i-(nr+1)):(d+(i-(nr+1)))+nb, i) = B;
    end
end

%% Solution of the Diophantine equation to find S, R polynomials 
theta = inv(M)*(V-Y);
R = [1; theta(1:nr)]';
S = [theta((nr+1):end)]';

% 
% R
% S