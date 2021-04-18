clear
% ��������� �����-��������������
Fs = 1000; % ������� �������������
Ts = 1/Fs; % ��� �� �������
% ����� ���������� �������
duration = [0 10]; % ��������� � �������� ����� �������

% ��������� �������� (������ �����������)
L = 1000; ro = 2700; d = 114;

% ���������� ������������ �-���, ����������� � W_line.m, 
[W_line, w, t, FftL] = W_line(Fs, ro, d, L, duration, 'P');

% ����� ������. ��� ��������� ��. � signal_type.m
Signal = signal_type(t, 'sine');

%--------------------------�����-��������������---------------------------
%-------------------------------------------------------------------------
% ����� ������� ��������� (��� ����)
original_fft = fft(Signal, FftL);

% ����� ��������� ������� ����� ����������� ����� (���������)
Fft_S_out = W_line .* (original_fft);

% ������� �� ��������� ������� � ������� ������������ �-��� W_line (���� ���)
% ��� ������� ������: �������������� �������� ������� �� ���������
fft_restored = Fft_S_out ./ W_line;


% ���������� �� ��������� ����� **(�� �����������, �������� ������������
% ������������� fft (��� ����, �������� ��� �� ������� ������)
restored_signal = ifft(fft_restored, FftL);


% �������� ������, ��������� ����� (����������� �� �����������) 
S_out = ifft(Fft_S_out);
%-------------------------------------------------------------------------
% ������ ��� (� ������ �������� ������ ������� �������� � ����)
err = rand(1,length(t));
err = err - mean(err);
err = err * real(max(S_out))*0.2 / std(err); % ���/������ ������
%-------------------------------------------------------------------------

% ��������� ��� � ��������� �������, ��������� �����
S_out_noise = S_out + err;

% ����������� ������� � ��������� �����
Fft_S_out_noise = fft(S_out_noise, FftL);

%-------------------------��������� ������--------------------------------
Fft_S_out_filtered = Fft_S_out_noise .* ...
    bandpass_filter(Fs, FftL, 1024, 2, 'yes');
%-------------------------------------------------------------------------

% ��������������� ������� ������ �� ������� ��������� 
% � ������� ������������ �������
% �� ���������� ��������� ������!
fft_restored_signal_raw = Fft_S_out_noise ./ W_line;
restored_signal_noise_raw = ifft(fft_restored_signal_raw, FftL);
restored_signal_noise_raw = real(restored_signal_noise_raw); % ����������� �� ������ ���������

% ��� ��� ��������������� ������� ������, �� ������ �������� � �����-������
% ��������� ������� ��������� ������ (��. ����, �-��� bandpass_filter.m)
fft_restored_signal_filtered = Fft_S_out_filtered ./ W_line;
restored_signal_filtered = ifft(fft_restored_signal_filtered, FftL);
restored_signal_filtered = real(restored_signal_filtered); % ����������� �� ������ ���������

%---------����������� !���������������� � �������������o! �������----------
%--------------------������ ����������-�����-------------------------------
Sav_Gol_filtered = sgolayfilt(restored_signal_filtered, 3, 251);

%-------------�������������� ������������ MA-������----------------------
MA_filtered = MA_filter(Sav_Gol_filtered, 50);


%------------------------������������ �����������-------------------------
% �������
F = Fs*(0:FftL/2-1)/FftL; % �����

% 1) ������ �������, ��� ������ �� ����
figure();
plot(t, Signal);
title('������� ������');
xlabel('t, s');
ylabel('A');

% 2) �������������� �������� ������� �� ��������� ��� ������� ����
% ��� ������������� ���������� �������
figure();
plot(t, restored_signal_noise_raw);
title('������� ������ �� ������� ���������, ��� �������');
xlabel('t, s');
ylabel('A');

% 3) �������������� �������� ������� �� ��������� ��� ������� ����
% � �������������� ���������� �������
figure();
plot(t, restored_signal_filtered);
title('������� ������ �� ������� ���������, � ��������');
xlabel('t, s');
ylabel('A');

% 4) ����� �������� ��������
figure();
% �������� ������, �������������� �� ����������� � ��������� ������ 
% (���� ���)
subplot(2,1,1);
plot(t, real(S_out));
title('������ �������� ������');
xlabel('t, s');
ylabel('A');

% �������� ������ ��� ������� ���� 
subplot(2,1,2);
plot(t, real(S_out_noise));
title('������ �������� ������');
xlabel('t, s');
ylabel('A');