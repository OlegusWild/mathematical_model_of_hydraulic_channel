% creating moving average filter
counts_per_frame = 250;
coef25MA = ones(1, counts_per_frame)/counts_per_frame;
MA_sig = filter(coef25MA, 1, cubic_test_filt);

figure('Name', '����������� ���������������� ��������� ������� MA', 'NumberTitle', 'off');
plot(t, MA_sig);
title('������ �������� ������ ����� ��������������� �����������');
xlabel('t, s');
ylabel('A');