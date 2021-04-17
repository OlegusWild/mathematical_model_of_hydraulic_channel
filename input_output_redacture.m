clear
% ��������� �����-��������������
Fs = 1000; % ������� �������������
Ts = 1/Fs; % ��� �� �������

% ����� ���������� �������
%duration = [0 30]; % ��������� � �������� ����� �������
duration = [0 10]; % ��������� � �������� ����� �������


FftL = 2^nextpow2(length(duration(1):Ts:duration(2))); % ���-�� ����� ��������������
t = linspace(duration(1),duration(2), FftL); % � ������������ � ���-��� ����� ������ ��������� �����
w = 0:Fs/FftL:Fs - Fs/FftL; % "��������" �������� ������� >= 0

% ��������� ������������ ������� ��������
    % �����������
    a = 0.4; % 1/�
    % �������
    L = 1000; % �
    % ������
    K = 4*10^(-8); % �^3/(��*�)
    K_l = 2.2 * 10^9; % ��
    % ������ ���� ��� ��������� �����
    E = 2 * 10^11; % ��
    % ���������� ������� ������������
    T = 0.01; % 1/�
    % �������������� �������������
     R = 1e6; % ��/(�^3/�)^2
      R_1 = R/2;
      R_L = R/2;
    % ��������� �������� ��������
    ro = 2700; % ��/�^3
    K_pr = K_l / (1 + a*K_l/E);
    c = sqrt(K_pr / ro); % �������� ����� � ��������, �/�
    % ������� �����
    d = 114; % ��
    % ������� ������� �����
    f = pi*d^2/4 *10^(-6); % �^2
    
    % ������ ������������ �-���
    v = 1/c*(w/2).^0.5 .* ((sqrt(w.^2+4*a^2)-w).^0.5 + 1i*(sqrt(w.^2+4*a^2)+w).^0.5);

    ro_g = c/f * (1 + 4*a^2 ./ w.^2).^0.25 .* exp(-1i*0.5*atan(2*a./w));

    W = (1i * w * K) ./ (1 + 1i*w.*T);
    

    % �� ��������
    W_line = ( (1+2*R*W) .* cosh(v*L) + (W*ro .* ro_g + ...
    2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);

    W_line(1) = 1; % ������ ��������� NaN �� ������� �������

    % �� �������
%     W_line = W .* ( (1+2*R*W) .* cosh(v*L) + (W*ro.*ro_g + ...
%     2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);
%     W_line = W_line * 10^8 * 0.8;

   

  


% --------------������� ������� ��������� ������--------------------------
% order = 64; % ��� �������� ���� 2048
% cut_off = (0.5) / (Fs/2); % cut_off fr-cy in pi units (0.5 Hz)
% cut_off_2 = (2.5) / (Fs/2);
% h = fir1(order, [cut_off cut_off_2]);

% ����� �������
% fft_h = fft(h, FftL);
% fft_h_A = abs(fft_h) * 2 ./ FftL;
% fft_h_A(1) = fft_h_A(1) / 2;


% --------------------------����������� ������----------------------------
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
% T_s = 2;
% tau = 1;
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

% ������� �� ��������� ������� � ������� ������������ �-��� W_line (���� ���)
fft_restored = Fft_S_out ./ W_line;
restored_signal = ifft(fft_restored, FftL);

% Output signal in time domain 
S_out = ifft(Fft_S_out);
%----------------------------
% ������ ���
err = rand(1,length(t));
err = err - mean(err);
err = err * real(max(S_out))*0.2 / std(err);
% ��������� ��� � ��������� �������
S_out_noise = S_out + err;
% ����������� ������� � ��������� �����
Fft_S_out_noise = fft(S_out_noise, FftL);


%-------------------------------filter apply------------------------------
% Fft_S_out_noise = Fft_S_out_noise .* fft_h;
%-------------------------------------------------------------------------


% ��������������� ������� ������ �� ������� ��������� 
% � ������� ������������ �������
fft_restored_signal = Fft_S_out_noise ./ W_line;
restored_signal_noise = ifft(fft_restored_signal, FftL);
restored_signal_noise = real(restored_signal_noise); % ����������� �� ������ ���������


%-------------------------������������ ������---------------------------
% cubic_test_filt = sgolayfilt(restored_signal_noise, 3, 251);

%-------------�������������� ������������ MA-������----------------------
% counts_per_frame = 50;
% coef25MA = ones(1, counts_per_frame)/counts_per_frame;
% MA_sig = filter(coef25MA, 1, restored_signal_noise);



% ��� ��� ������ ��������� �������
original_fft_A = abs(original_fft);
original_fft_A = 2*original_fft_A ./ FftL;
original_fft_A(1) = original_fft_A(1)/2;

% ��� ��� ������ ��������� �������
Fft_out_A=abs(Fft_S_out);% ��������� �������������� ����� �������
Fft_out_A=2*Fft_out_A./FftL;% ���������� ������� �� ���������
Fft_out_A(1)=Fft_out_A(1)/2;% ���������� ���������� ������������ � �������

% �������
F = Fs*(0:FftL/2-1)/FftL; % shifting


%------����������� �������� � ����������� ��������� �����������----------
% sig_diff = abs( real(restored_signal_noise) - Signal);
% sig_diff_squares = sig_diff .^2;
% sq_err = sqrt(sum(sig_diff_squares));

%-----------------------filter ���-----------------------------------------
% fft_h_A = abs(fft_h);
% fft_h_A = 2*fft_h_A ./ FftL;
% fft_h_A(1) = fft_h_A(1)/2;


%-----------------------------------���--------------------------------
% figure('Name', '��� ��� �������� � ��������� ��������', 'NumberTitle', 'off');
% subplot(2,1,1);
% plot(F,original_fft_A(1:length(F)));
% title('��� �������� �������');
% xlabel('F, Hz');
% ylabel('A');
% 
% subplot(2,1,2);
% plot(F, Fft_out_A(1:length(F)));
% title('��� ��������� �������');
% xlabel('F, Hz');
% ylabel('A');

%-----------------------------------��� �������------------------------
% figure('Name', '��� ��������� �������')
% plot(F, fft_h_A(1:length(F)));
% title('��� ��������� �������');
% xlabel('F, Hz');
% ylabel('A');


% ������ �������, ��� ������ �� ���� ���
figure('Name', '��������� ������� ������', 'NumberTitle', 'off');
% check for original signal
S_orig = ifft(original_fft);
plot(t,S_orig);
title('������� ������');
xlabel('t, s');
ylabel('A');


% ----------------getting signals in time domain---------------------
figure('Name', '������� �������, ��������������� �� ������� � ������� ���������', 'NumberTitle', 'off');
subplot(2,1,1);
% ������� ������, ��������������� �� ������� ��������� � ������� W_line
plot(t, real(restored_signal));
title('������� ������ �� ������� ���������');
xlabel('t, s');
ylabel('A');
% ylim([-0.2 1.2]);

% ������� ������, ��������������� �� ��������� ������� � ������� W_line
subplot(2,1,2);
plot(t, (restored_signal_noise));
title('������� ������ �� ������� ���������');
xlabel('t, s');
ylabel('A');

%-----------------------����� �������� "������������" ��������-----------
figure('Name', '�������� �������: ������ � � �����', 'NumberTitle', 'off');
% �������� "������������" ������ �� �����������
subplot(2,1,1);
plot(t, real(S_out));
title('������ �������� ������');
xlabel('t, s');
ylabel('A');

% ��������, �� ��� � � �����
subplot(2,1,2);
plot(t, real(S_out_noise));
title('������ �������� ������');
xlabel('t, s');
ylabel('A');


% %-----------------------������������������ ������-----------------------
% figure('Name', '������������������ ������', 'NumberTitle', 'off');
% plot(A_err, sq_err);
% title('������ �������� ������');
% xlabel('t, s');
% ylabel('A');

%--------------------���������� ������������ ��������----------------------
% figure('Name', '����������� ���������������� ��������� �������', 'NumberTitle', 'off');
% plot(t, cubic_test_filt);
% title('������ �������� ������ ����� ��������������� �����������');
% xlabel('t, s');
% ylabel('A');
% 
% %------------------�������������� ����������� MA-��������-----------------
figure('Name', '����������� ���������������� ��������� ������� ��� ������ MA-�������', 'NumberTitle', 'off');
plot(t, MA_sig);
title('������ �������� ������ ����� ��������������� ����������� ��� ������ MA-�������');
xlabel('t, s');
ylabel('A');