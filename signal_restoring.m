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
% ������ ��� (� ������ �������� ������ ������� �������� � ����) (���������)
err = rand(1,length(t));
err = err - mean(err);
err1 = err * real(max(S_out))*0.01 / std(err); % ���/������ ������
err2 = err * real(max(S_out))*0.03 / std(err);
err3 = err * real(max(S_out))*0.1 / std(err);
%-------------------------------------------------------------------------

% ��������� ��� � ��������� �������, ��������� ����� (���������)
S_out_noise1 = S_out + err1;
S_out_noise2 = S_out + err2;
S_out_noise3 = S_out + err3;

% ����������� ������� � ��������� ����� (���������)
Fft_S_out_noise1 = fft(S_out_noise1, FftL);
Fft_S_out_noise2 = fft(S_out_noise2, FftL);
Fft_S_out_noise3 = fft(S_out_noise3, FftL);

%-------------------------��������� ������--------------------------------
%(���������)
Fft_S_out_filtered1 = Fft_S_out_noise1 .* ...
    bandpass_filter(Fs, FftL, 512, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 512, 5, 'no')));
    % ����� �������, ��� ��� ���������� ������� �� �������!
Fft_S_out_filtered2 = Fft_S_out_noise2 .* ...
    bandpass_filter(Fs, FftL, 512, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 512, 5, 'no')));
Fft_S_out_filtered3 = Fft_S_out_noise3 .* ...
    bandpass_filter(Fs, FftL, 512, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 512, 5, 'no')));
%-------------------------------------------------------------------------

% ��������������� ������� ������ �� ������� ��������� (���������)
% � ������� ������������ �������
% �� ���������� ��������� ������!
fft_restored_signal_raw1 = Fft_S_out_noise1 ./ W_line;
restored_signal_noise_raw1 = ifft(fft_restored_signal_raw1, FftL);
restored_signal_noise_raw1 = real(restored_signal_noise_raw1); % ����������� �� ������ ���������

fft_restored_signal_raw2 = Fft_S_out_noise2 ./ W_line;
restored_signal_noise_raw2 = ifft(fft_restored_signal_raw2, FftL);
restored_signal_noise_raw2 = real(restored_signal_noise_raw2);

fft_restored_signal_raw3 = Fft_S_out_noise3 ./ W_line;
restored_signal_noise_raw3 = ifft(fft_restored_signal_raw3, FftL);
restored_signal_noise_raw3 = real(restored_signal_noise_raw3);

% ��� ��� ��������������� ������� ������, �� ������ �������� � �����-������
% (���������)
% ��������� ������� ��������� ������ (��. ����, �-��� bandpass_filter.m)
fft_restored_signal_filtered1 = Fft_S_out_filtered1 ./ W_line;
restored_signal_filtered1 = ifft(fft_restored_signal_filtered1, FftL);
restored_signal_filtered1 = real(restored_signal_filtered1); % ����������� �� ������ ���������

fft_restored_signal_filtered2 = Fft_S_out_filtered2 ./ W_line;
restored_signal_filtered2 = ifft(fft_restored_signal_filtered2, FftL);
restored_signal_filtered2 = real(restored_signal_filtered2);

fft_restored_signal_filtered3 = Fft_S_out_filtered3 ./ W_line;
restored_signal_filtered3 = ifft(fft_restored_signal_filtered3, FftL);
restored_signal_filtered3 = real(restored_signal_filtered3);

%---------����������� !���������������� � �������������o! �������----------
%--------------------������ ����������-�����-------------------------------
% Sav_Gol_filtered1 = sgolayfilt(restored_signal_filtered1, 3, 251);
% Sav_Gol_filtered2 = sgolayfilt(restored_signal_filtered2, 3, 251);
% Sav_Gol_filtered3 = sgolayfilt(restored_signal_filtered3, 3, 251);

%-------------�������������� ������������ MA-������----------------------
%MA_filtered = MA_filter(Sav_Gol_filtered, 50);

%--------------------------����� ������ ��� �������-----------------------
figure();
subplot(3,2,1);
plot(t, restored_signal_noise_raw1);
title('������� ������ �� ������� ���������, ��� �������');
xlabel('t, s');
ylabel('A');

subplot(3,2,3);
plot(t, restored_signal_noise_raw2);
xlabel('t, s');
ylabel('A');

subplot(3,2,5);
plot(t, restored_signal_noise_raw3);
xlabel('t, s');
ylabel('A');

subplot(3,2,2);
plot(t, restored_signal_filtered1);
title('������� ������ �� ������� ���������, � ��������');
xlabel('t, s');
ylabel('A');

subplot(3,2,4);
plot(t, restored_signal_filtered2);
xlabel('t, s');
ylabel('A');

subplot(3,2,6);
plot(t, restored_signal_filtered3);
xlabel('t, s');
ylabel('A');

% %------------------------������������ �����������-------------------------
% % �������
% F = Fs*(0:FftL/2-1)/FftL; % �����
% 
% 1) ������ �������, ��� ������ �� ����
figure();
plot(t, Signal);
title('������� ������');
xlabel('t, s');
ylabel('A');
% 
% % 2) �������������� �������� ������� �� ��������� ��� ������� ����
% % ��� ������������� ���������� �������
% figure();
% plot(t, restored_signal_noise_raw);
% title('������� ������ �� ������� ���������, ��� �������');
% xlabel('t, s');
% ylabel('A');
% 
% % 3) �������������� �������� ������� �� ��������� ��� ������� ����
% % � �������������� ���������� �������
% figure();
% plot(t, restored_signal_filtered);
% title('������� ������ �� ������� ���������, � ��������');
% xlabel('t, s');
% ylabel('A');

% % 4) ����� �������� ��������
% figure();
% % �������� ������, �������������� �� ����������� � ��������� ������ 
% % (���� ���)
% subplot(2,1,1);
% plot(t, real(S_out));
% title('������ �������� ������');
% xlabel('t, s');
% ylabel('A');
% 
% % �������� ������ ��� ������� ���� 
% subplot(2,1,2);
% plot(t, real(S_out_noise));
% title('������ �������� ������');
% xlabel('t, s');
% ylabel('A');