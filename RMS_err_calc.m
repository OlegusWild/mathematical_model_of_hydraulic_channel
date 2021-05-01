clear
% Parameters of Fft
Fs = 1000; % sampling, depends on apparatus
Ts = 1/Fs; % timestep

%duration = [0 10]; % initial and final points of time domain
duration = [0 10]; % initial and final points of time domain


% ��������� �������� (������ �����������)
L = 1000; ro = 2700; d = 114;

% ���������� ������������ �-���, ����������� � W_line.m, 
[W_line, w, t, FftL] = W_line(Fs, ro, d, L, duration, 'P');

% ����� ������. ��� ��������� ��. � signal_type.m
Signal = signal_type(t, 'harmonical');

%--------------------------�����-��������������---------------------------
%-------------------------------------------------------------------------
% ����� ������� ��������� (��� ����)
original_fft = fft(Signal, FftL);

% ����� ��������� ������� ����� ����������� ����� (���������)
Fft_S_out = W_line .* (original_fft);

% �������� ������, ��������� ����� (����������� �� �����������) 
S_out = ifft(Fft_S_out);
k = 1;
% ��������� ��������� ����, ������� ����������� � ��������� ������� (����� ���������������)

    for j=0:0.001:0.05
        % ������ ���
        err = rand(1,length(t));
        err = err - mean(err);
        err = err * real(max(S_out)) *j / std(err);
        % ��������� ��� � ��������� �������
        S_out_noise = S_out + err;
        % ����������� ������� � ��������� �����
        Fft_S_out_noise = fft(S_out_noise, FftL);

%-------------------------------------���� ��� ����������-----------------
        % ��������������� ������� ������ �� ������� ��������� 
        % � ������� ������������ �������
%         fft_restored_signal = Fft_S_out_noise ./ W_line;
%         restored_signal_noise = ifft(fft_restored_signal, FftL);



%---------------------------------���� ����, ���� � �����������-----------

        %-------------------------------filter apply------------------------------
        Fft_S_out_filtered = Fft_S_out_noise .* ...
    bandpass_filter(Fs, FftL, 64, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 64, 5, 'no')));
    % ����� �������, ��� ��� ���������� ������� �� �������!
        %-----------------------------------------------------------------
        % ��� ��� ��������������� ������� ������, �� ������ �������� � �����-������
        % ��������� ������� ��������� ������ (��. ����, �-��� bandpass_filter.m)
        fft_restored_signal_filtered = Fft_S_out_filtered ./ W_line;
        restored_signal_filtered = ifft(fft_restored_signal_filtered, FftL);
        restored_signal_filtered = real(restored_signal_filtered); % ����������� �� ������ ���������


        
        
        %------����������� �������� � ����������� ��������� �����������----------
        sig_diff = abs( (restored_signal_filtered) - Signal);
%         sig_diff = abs( real(restored_signal_noise) - Signal);
        sig_diff_squares = sqrt(sig_diff .^2);
        %ind = int8(j*100);

        sq_err(k) = sum(sig_diff_squares) / length(Signal); % ����� �� ���-�� ��������
        k = k+1;
    end


    % � ����� ����� ����� -
    %-----------------------������������������ ������-----------------------
    
x = (0:0.001:0.05)*100;
    
figure('Name', '������������� ������������������ ������', 'NumberTitle', 'off');

hold on;
plot(x, sq_err/(max(Signal)));
plot((0:0.005:0.05)*100, sq_err(1:5:end)/(max(Signal)), 'k.', 'MarkerSize', 15)

title('������������� ������������������ ������');
xlabel('���/������, %');
%ylabel('������������� ��.��.������, 1/A^2'); % 10^5*