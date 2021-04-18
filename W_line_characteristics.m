clear
% Parameters of Fft
Fs = 1000; % sampling, depends on apparatus
Ts = 1/Fs; % timestep

%duration = [0 10]; % initial and final points of time domain
duration = [0 10]; % initial and final points of time domain


FftL = 2^nextpow2(length(duration(1):Ts:duration(2))); % number of points fot FFT
t = linspace(duration(1),duration(2), FftL);
w = 0:Fs/FftL:Fs - Fs/FftL;

% Частоты
F = Fs*(0:FftL/2-1)/FftL; % shifting


% figure('Name', 'АЧХ передаточной функции', 'NumberTitle', 'off');
for i=2:7
hold on;
% parameters for transfer-function for the drill line
    a = 0.4; % 1/s
    L = 1000*i; % m
    K = 4*10^(-8); % m^3/Pa/s
    K_l = 2.2 * 10^9; % Pa
    E = 2 * 10^11; % Pa
    T = 0.01; % 1/s
     R = 1e6; % ?
      R_1 = R/2;
      R_L = R/2;
    
    ro = 2700; % kg/m^3
    K_pr = K_l / (1 + a*K_l/E);
    c = sqrt(K_pr / ro);

    d = 114; % mm
    f = pi*d^2/4 *10^(-6); % m^2

    
    % matrix values
    v = 1/c*(w/2).^0.5 .* ((sqrt(w.^2+4*a^2)-w).^0.5 + 1i*(sqrt(w.^2+4*a^2)+w).^0.5);

    ro_g = c/f * (1 + 4*a^2 ./ w.^2).^0.25 .* exp(-1i*0.5*atan(2*a./w));

    W = (1i * w * K) ./ (1 + 1i*w.*T);
    

    % pressure
    W_line = ( (1+2*R*W) .* cosh(v*L) + (W*ro .* ro_g + ...
    2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);

    W_line(1) = 1; % тут NaN!

    % Q
%     W_line = W .* ( (1+2*R*W) .* cosh(v*L) + (W*ro.*ro_g + ...
%     2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);
%     W_line = W_line * 10^8 * 0.8;

     
    % Строим АЧХ
 plot(F, abs(W_line(1:length(F))),'--', 'linewidth', 1.2);
    title('АЧХ передаточной функции по давлению');
    xlabel('F, Гц');
    ylabel('A');   
    xlim([0 3])
%    
    % Для построения ФЧХ
%   plot(F, angle(W_line(1:length(F)))/pi,'--', 'linewidth', 1);
%     title('ФЧХ передаточной функции по расходу');
%     xlabel('F, Гц');
%     ylabel('ед. \pi');   
%     xlim([0 3])
%     
    
    




    
end


% legend(string(2e3)+' м', string(3e3)+' м', string(4e3)+' м', string(5e3)+' м', string(6e3)+' м', string(7e3)+' м')