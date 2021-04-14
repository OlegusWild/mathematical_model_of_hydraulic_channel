clear
% Parameters of Fft
Fs = 1000; % sampling, depends on apparatus
Ts = 1/Fs; % timestep

%duration = [0 10]; % initial and final points of time domain
duration = [0 10]; % initial and final points of time domain


FftL = 2^nextpow2(length(duration(1):Ts:duration(2))); % number of points for FFT
t = linspace(duration(1),duration(2), FftL);
w = 0:Fs/FftL:Fs - Fs/FftL;

% parameters for transfer-function for the drill line
    a = 0.4; % 1/s
    L = 1000; % m
    K = 4*10^(-8); % m^3/Pa/s
    K_l = 2.2 * 10^9; % Pa
    E = 2 * 10^11; % Pa
    T = 0.01; % 1/s
     R = 1e10; % ?
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

    W_line(1) = -0.006; % ��� NaN!

    % Q
%     W_line = W .* ( (1+2*R*W) .* cosh(v*L) + (W*ro.*ro_g + ...
%     2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);
%     W_line = W_line * 10^8 * 0.8;

    % abs_W_line = sqrt((real(W_line)).^2 + (imag(W_line)).^2);


% �� ������, ���� ����� ����� ��������� ��� �������������� �������
% --------------������� ������� �������� ������--------------------------
% order = 32;
% cut_off = (0.5) / (Fs/2); % cut_off fr-cy in pi units (0.5 Hz)
% cut_off_2 = (7) / (Fs/2);
% h = fir1(order, [cut_off cut_off_2]);
% % ����� �������
% fft_h = fft(h, FftL);
% fft_h_A = abs(fft_h) * 2 ./ FftL;
% fft_h_A(1) = fft_h_A(1) / 2;


% --------------------------experimental signal--------------------------
% ����� ������� ��� �������
Signal = sin(2*pi*1.7*t);
% Signal = 5*sin(2*pi*1.7*t) + 5*sin(2*pi*10*t) + 2*sin(2*pi*0.7*t);

% ? how ? 5) ������� ������� ��������� 0 1 0 �� ������� 30 ���
% Signal = [];
% Signal(1:length(t)/3) = 7*sin(2*pi*t(1:length(t)/3));
% Signal(length(t)/3+1:2*length(t)/3) = 7*sin(2*pi*t(length(t)/3+1:2*length(t)/3) - pi/4);
% Signal(2*length(t)/3 + 1:length(t)) = 7*sin(2*pi*t(2*length(t)/3 + 1:length(t)));

% 6) ������� ��������� ��������� 0 1 0
% Signal = [];
% Signal(1:length(t)/3) = 7*sin(2*pi*t(1:length(t)/3));
% Signal(length(t)/3+1:2*length(t)/3) = 5*sin(2*pi*4*t(length(t)/3+1:2*length(t)/3));
% Signal(2*length(t)/3 + 1:length(t)) = 7*sin(2*pi*t(2*length(t)/3 + 1:length(t)));


% 1) ������� ������������� ��������
% T_s = 1;
% tau = 0.2;
% Signal = rem(t,T_s)<tau; % �������������� �������

% 2) �������� � "1", �������� � "0"
% Signal(1:length(t)/2) = rem(t(1:length(t)/2),T_s)<tau; % �������������� �������
% Signal(length(t)/2+1:length(t)) = zeros([1 length(t)/2]);

% 3) ��� �� ��������
% tau = 0.2;
% here we make a pulse of tau duration
% [row, col] = find(t<5+tau/2 & t>5-tau/2);
% Signal(col(1):col(end)) = 1; % tooth
% Signal(1:col(1)-1) = 0;
% Signal(col(end)+1:length(t)) = 0;

% 4) SyncSequence_1
% tau = 2;
% T = 2*tau;
% duration 30 s
%������� ���������� ���������������� ������������������ � �����
% [row, col] = find(t>=3+tau/2 & t<=23+tau/2);
% Signal(col(1):col(end)) = (rem(t(col(1):col(end)),T)<tau);
% Signal(1:col(1)-1) = 0;
% % ����� �������� ���������
% [row,col] = find(t>=23+tau/2 & t<=26+tau/2);
% Signal(col(1):col(end)) = (rem(t(col(1):col(end)),T/2)<tau/2);
% Signal(col(end)+1:length(t)) = 0;

% ����� ������� ��������� (��� ����)
original_fft = fft(Signal, FftL);

% ����� ��������� ������� ����� ����������� �����
Fft_S_out = W_line .* (original_fft);

% Output signal in time domain 
S_out = ifft(Fft_S_out);
%----------------------------
k = 1;
% ��������� ��������� ����

    for j=0:0.001:0.05
        % ������ ���
        err = rand(1,length(t));
        err = err - mean(err);
        err = err * real(max(S_out)) *j / std(err);
        % ��������� ��� � ��������� �������
        S_out_noise = S_out + err;
        % ����������� ������� � ��������� �����
        Fft_S_out_noise = fft(S_out_noise, FftL);

        % ��������������� ������� ������ �� ������� ��������� 
        % � ������� ������������ �������
        fft_restored_signal = Fft_S_out_noise ./ W_line;
        restored_signal_noise = ifft(fft_restored_signal, FftL);

        %-------------------------------filter apply------------------------------
        % Fft_S_out = Fft_S_out .* fft_h;
        %-------------------------------------------------------------------------

        %------����������� �������� � ����������� ��������� �����������----------
        sig_diff = abs( real(restored_signal_noise) - Signal);
        sig_diff_squares = sqrt(sig_diff .^2);
        %ind = int8(j*100);

        sq_err(k) = sum(sig_diff_squares) / length(Signal); % ����� �� ���-�� ��������
        k = k+1;
    end


    % � ����� ����� ����� -
    %-----------------------������������������ ������-----------------------
    

    
figure('Name', '������������� ������������������ ������', 'NumberTitle', 'off');

hold on;
plot((0:0.001:0.05)*100, sq_err/(max(Signal))); % normalized!!! /1e5
plot((0:0.005:0.05)*100, sq_err(1:5:end)/(max(Signal)), 'k.', 'MarkerSize', 15)

title('������������� ������������������ ������');
xlabel('���/������, %');
%ylabel('������������� ��.��.������, 1/A^2'); % 10^5*