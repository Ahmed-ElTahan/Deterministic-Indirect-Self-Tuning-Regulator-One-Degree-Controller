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
