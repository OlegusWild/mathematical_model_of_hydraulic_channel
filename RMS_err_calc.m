clear
% Parameters of Fft
Fs = 1000; % sampling, depends on apparatus
Ts = 1/Fs; % timestep

%duration = [0 10]; % initial and final points of time domain
duration = [0 10]; % initial and final points of time domain


% Параметры скважины (только варьируемые)
L = 1000; ro = 2700; d = 114;

% Используем передаточную ф-цию, формируемую в W_line.m, 
[W_line, w, t, FftL] = W_line(Fs, ro, d, L, duration, 'P');

% Задаём сигнал. Все доступные см. в signal_type.m
Signal = signal_type(t, 'harmonical');

%--------------------------Фурье-преобразования---------------------------
%-------------------------------------------------------------------------
% Образ сигнала идеальный (без шума)
original_fft = fft(Signal, FftL);

% Образ выходного сигнала после прохождения линии (симуляция)
Fft_S_out = W_line .* (original_fft);

% Выходной сигнал, временной домен (регистрация на поверхности) 
S_out = ifft(Fft_S_out);
k = 1;
% Варьируем амплитуду шума, который добавляется к выходному сигналу (перед восстановлением)

    for j=0:0.001:0.05
        % Создаём шум
        err = rand(1,length(t));
        err = err - mean(err);
        err = err * real(max(S_out)) *j / std(err);
        % Добавляем шум к выходному сигналу
        S_out_noise = S_out + err;
        % Преобразуем обратно в частотный домен
        Fft_S_out_noise = fft(S_out_noise, FftL);

%-------------------------------------блок без фильтрации-----------------
        % Восстанавливаем входной сигнал по шумному выходному 
        % с помощью передаточной функции
%         fft_restored_signal = Fft_S_out_noise ./ W_line;
%         restored_signal_noise = ifft(fft_restored_signal, FftL);



%---------------------------------этот блок, если с фильтрацией-----------

        %-------------------------------filter apply------------------------------
        Fft_S_out_filtered = Fft_S_out_noise .* ...
    bandpass_filter(Fs, FftL, 64, 5, 'no') .* exp(-1i*w.*angle(bandpass_filter(Fs, FftL, 64, 5, 'no')));
    % Важно помнить, что ФЧХ полосового фильтра не нулевое!
        %-----------------------------------------------------------------
        % Ещё раз восстанавливаем входной сигнал, но теперь применяя к Фурье-образу
        % выходного сигнала полосовой фильтр (см. выше, ф-ция bandpass_filter.m)
        fft_restored_signal_filtered = Fft_S_out_filtered ./ W_line;
        restored_signal_filtered = ifft(fft_restored_signal_filtered, FftL);
        restored_signal_filtered = real(restored_signal_filtered); % Избавляемся от мнимых дребезгов


        
        
        %------расхождение сигналов с увеличением амплитуды погрешности----------
        sig_diff = abs( (restored_signal_filtered) - Signal);
%         sig_diff = abs( real(restored_signal_noise) - Signal);
        sig_diff_squares = sqrt(sig_diff .^2);
        %ind = int8(j*100);

        sq_err(k) = sum(sig_diff_squares) / length(Signal); % делим на кол-во отсчётов
        k = k+1;
    end


    % В итоге после цикла -
    %-----------------------среднеквадратичная ошибка-----------------------
    
x = (0:0.001:0.05)*100;
    
figure('Name', 'Относительная среднеквадратичная ошибка', 'NumberTitle', 'off');

hold on;
plot(x, sq_err/(max(Signal)));
plot((0:0.005:0.05)*100, sq_err(1:5:end)/(max(Signal)), 'k.', 'MarkerSize', 15)

title('Относительная среднеквадратичная ошибка');
xlabel('шум/сигнал, %');
%ylabel('Относительная Ср.кв.ошибка, 1/A^2'); % 10^5*