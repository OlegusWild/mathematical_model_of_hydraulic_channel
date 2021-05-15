clear
% Параметры Фурье-преобразования
Fs = 1000; % частота дискретизации
Ts = 1/Fs; % шаг по времени
% Выбор временного отрезка
duration = [0 10]; % начальная и конечная точки отрезка

% Параметры скважины (только варьируемые)
L = 1000; ro = 2700; d = 114;

% Используем передаточную ф-цию, формируемую в W_line.m, 
[W_line, w, t, FftL] = W_line(Fs, ro, d, L, duration, 'P');

% Задаём сигнал. Все доступные см. в signal_type.m
Signal = signal_type(t, 'sine');

%--------------------------Фурье-преобразования---------------------------
%-------------------------------------------------------------------------
% Образ сигнала идеальный (без шума)
original_fft = fft(Signal, FftL);

% Образ выходного сигнала после прохождения линии (симуляция)
Fft_S_out = W_line .* (original_fft);

% Входной по выходному сигналу с помощью передаточной ф-ции W_line (шума нет)
% Это братная задача: восстановление входного сигнала по выходному
fft_restored = Fft_S_out ./ W_line;


% Возвращаем во временной домен **(не обязательно, проверка правильности
% использования fft (нет шума, получаем тот же входной сигнал)
restored_signal = ifft(fft_restored, FftL);


% Выходной сигнал, временной домен (регистрация на поверхности) 
S_out = ifft(Fft_S_out);
%-------------------------------------------------------------------------
% Создаём шум (в датчик давления помимо сигнала попадают и шумы) (почистить)
err = rand(1,length(t));
err = err - mean(err);
err1 = err * real(max(S_out))*0.01 / std(err); % шум/сигнал опцион
err2 = err * real(max(S_out))*0.03 / std(err);
err3 = err * real(max(S_out))*0.1 / std(err);
%-------------------------------------------------------------------------

% Добавляем шум к выходному сигналу, временной домен (почистить)
S_out_noise1 = S_out + err1;
S_out_noise2 = S_out + err2;
S_out_noise3 = S_out + err3;

% Преобразуем обратно в частотный домен (почистить)
Fft_S_out_noise1 = fft(S_out_noise1, FftL);
Fft_S_out_noise2 = fft(S_out_noise2, FftL);
Fft_S_out_noise3 = fft(S_out_noise3, FftL);

%-------------------------полосовой фильтр--------------------------------
%(почистить)
Fft_S_out_filtered1 = Fft_S_out_noise1 .* ...
    bandpass_filter(Fs, FftL, 512, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 512, 5, 'no')));
    % Важно помнить, что ФЧХ полосового фильтра не нулевое!
Fft_S_out_filtered2 = Fft_S_out_noise2 .* ...
    bandpass_filter(Fs, FftL, 512, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 512, 5, 'no')));
Fft_S_out_filtered3 = Fft_S_out_noise3 .* ...
    bandpass_filter(Fs, FftL, 512, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 512, 5, 'no')));
%-------------------------------------------------------------------------

% Восстанавливаем входной сигнал по шумному выходному (почистить)
% с помощью передаточной функции
% НЕ используем полосовой фильтр!
fft_restored_signal_raw1 = Fft_S_out_noise1 ./ W_line;
restored_signal_noise_raw1 = ifft(fft_restored_signal_raw1, FftL);
restored_signal_noise_raw1 = real(restored_signal_noise_raw1); % Избавляемся от мнимых дребезгов

fft_restored_signal_raw2 = Fft_S_out_noise2 ./ W_line;
restored_signal_noise_raw2 = ifft(fft_restored_signal_raw2, FftL);
restored_signal_noise_raw2 = real(restored_signal_noise_raw2);

fft_restored_signal_raw3 = Fft_S_out_noise3 ./ W_line;
restored_signal_noise_raw3 = ifft(fft_restored_signal_raw3, FftL);
restored_signal_noise_raw3 = real(restored_signal_noise_raw3);

% Ещё раз восстанавливаем входной сигнал, но теперь применяя к Фурье-образу
% (почистить)
% выходного сигнала полосовой фильтр (см. выше, ф-ция bandpass_filter.m)
fft_restored_signal_filtered1 = Fft_S_out_filtered1 ./ W_line;
restored_signal_filtered1 = ifft(fft_restored_signal_filtered1, FftL);
restored_signal_filtered1 = real(restored_signal_filtered1); % Избавляемся от мнимых дребезгов

fft_restored_signal_filtered2 = Fft_S_out_filtered2 ./ W_line;
restored_signal_filtered2 = ifft(fft_restored_signal_filtered2, FftL);
restored_signal_filtered2 = real(restored_signal_filtered2);

fft_restored_signal_filtered3 = Fft_S_out_filtered3 ./ W_line;
restored_signal_filtered3 = ifft(fft_restored_signal_filtered3, FftL);
restored_signal_filtered3 = real(restored_signal_filtered3);

%---------сглаживание !восстановленного и фильтрованногo! сигнала----------
%--------------------фильтр Савитского-Голея-------------------------------
% Sav_Gol_filtered1 = sgolayfilt(restored_signal_filtered1, 3, 251);
% Sav_Gol_filtered2 = sgolayfilt(restored_signal_filtered2, 3, 251);
% Sav_Gol_filtered3 = sgolayfilt(restored_signal_filtered3, 3, 251);

%-------------Дополнительный сглаживающий MA-фильтр----------------------
%MA_filtered = MA_filter(Sav_Gol_filtered, 50);

%--------------------------общий график для диплома-----------------------
figure();
subplot(3,2,1);
plot(t, restored_signal_noise_raw1);
title('Входной сигнал из шумного выходного, без фильтра');
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
title('Входной сигнал из шумного выходного, с фильтром');
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

% %------------------------ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ-------------------------
% % Частоты
% F = Fs*(0:FftL/2-1)/FftL; % сдвиг
% 
% 1) Просто смотрим, что подали на вход
figure();
plot(t, Signal);
title('Входной сигнал');
xlabel('t, s');
ylabel('A');
% 
% % 2) Восстановление входного сигнала из выходного при наличии шума
% % БЕЗ использования полосового фильтра
% figure();
% plot(t, restored_signal_noise_raw);
% title('Входной сигнал из шумного выходного, без фильтра');
% xlabel('t, s');
% ylabel('A');
% 
% % 3) Восстановление входного сигнала из выходного при наличии шума
% % С использованием полосового фильтра
% figure();
% plot(t, restored_signal_filtered);
% title('Входной сигнал из шумного выходного, с фильтром');
% xlabel('t, s');
% ylabel('A');

% % 4) Обзор выходных сигналов
% figure();
% % Выходной сигнал, регистрируемый на поверхности в идеальном случае 
% % (шума нет)
% subplot(2,1,1);
% plot(t, real(S_out));
% title('Чистый выходной сигнал');
% xlabel('t, s');
% ylabel('A');
% 
% % Выходной сигнал при наличии шума 
% subplot(2,1,2);
% plot(t, real(S_out_noise));
% title('Шумный выходной сигнал');
% xlabel('t, s');
% ylabel('A');