% �������� ���������� ������� �� �������� ������.
% �-��� ���������� �����-����� ������� � ��������� �����������:
% Fs - ������� �������������, FftL - ���-�� ����� �����-��������������
% order - ������� (��� ������������� ��������� � ������������� ���� ~2048)
% cut_off_freequencies - ������� �������; ���� ������� 1 ������� - 
% low-pass filter
% ����� plotting: 'yes' - �������� ������ ��� �������, 'no' - ���.

% ������ �������������:
% bandpass_filter(1000, 16384, 128, [5 100], 'yes');

function fft_h = bandpass_filter(Fs, FftL, order, cut_off_freequencies, plotting)
    if length(cut_off_freequencies) == 1
        cut_off = cut_off_freequencies / (Fs/2);
        h = fir1(order, cut_off); 
    end
    if length(cut_off_freequencies) == 2
        cut_off_1 = cut_off_freequencies(1) / (Fs/2);
        cut_off_2 = cut_off_freequencies(2) / (Fs/2);
        h = fir1(order, [cut_off_1 cut_off_2]); 
    end

    % ����� �������
    fft_h = fft(h, FftL);
    
    if plotting == "yes"
        fft_h_A = abs(fft_h) * 2 ./ FftL;
        fft_h_A(1) = fft_h_A(1) / 2;
        F = Fs*(0:FftL/2-1)/FftL;
        plot(F, fft_h_A(1:length(F)));
    end