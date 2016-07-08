clc; clear all; close all;

%% Inputs 

% input initialization
    Ts = 0.1; % sample time
    t = 0:Ts:30;
    
% Input 1 (damped sine wave)
%     b = 0.1;
%     a = 2;
%     w = 2*pi/2.5; %rad/s
%     uc = a*exp(-b*t).*sin(w*t);
     
%Input 2 (square wave)
    uc = square(t, 50)';
    
% Input 3 (step)
%     uc = ones(length(t), 1);
    
% Input 4 (ramdom)
%     uc = wgn(length(t), 1, 0);
    
% Input plot over time
    figure(1)
    plot(t, uc, 'LineWidth', 2);
    xlabel('Time (sec.)');
    ylabel('Input');
    title('Time - Input plot')
    grid on
    print(1,'-dpng','-r500','Input')
%% Transfer function preparation

    s = tf('s');
    
% System
    Gs = (s+1.5)*exp(-Ts*s)/(s-1)/(s+3);
   
% System z-transform and its properties
    Gz = c2d(Gs, Ts, 'zoh') % z-transformation
    [num, den] = tfdata(Gz); % get numerator and denomenator of the pulse tf
    num = cell2mat(num);
    den = cell2mat(den);
    na = length(den) -1;
    iodel = get(Gz, 'iodelay'); % input output delay
    d = max(find(num ==0))+iodel; % delay found in B matrix
    nb = length(num) - max(find(num ==0)) - 1;
    num = num(max(find(num ==0))+1:end); % numerator after removal of first zeros of delay
    nu = na + nb + 1; % number of unknowns

 % systems numerator and denomenator
    Asys = 1;
    for i = 1: na
        Asys(i+1) = den(i+1);
    end
        
    for i = na+1:nu
        Bsys(i-na) = num(i-na);
    end
           



%% Indirect Self-Tuning Regulator Fuction

% required closed loop polynomial
    % system
        z1 = exp(-Ts*2);
        z2 = exp(-Ts*4);
        Am = [1 -(z1+z2) z1*z2 ];

% Observer poles        
    % system
        A0 = [1 -0.4 0.04];

% Rrquired Characteristic Polynomial
    alpha = conv(Am, A0); 
    
% Indirect Self-Tuning Regulator (ISTR)
    [ theta, Gz_estm, Gc, Gcl] = ISTR( uc, Asys, Bsys, d, Ts, Am, A0) 

%% -------------------------------------------------------------------------------------------------------