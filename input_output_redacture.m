clear
% Параметры Фурье-преобразования
Fs = 1000; % частота дискретизации
Ts = 1/Fs; % шаг по времени

% Выбор временного отрезка
%duration = [0 30]; % начальная и конечная точки отрезка
duration = [0 10]; % начальная и конечная точки отрезка


FftL = 2^nextpow2(length(duration(1):Ts:duration(2))); % кол-во точек преобразования
t = linspace(duration(1),duration(2), FftL); % в соответствии с кол-вом точек строим временную шкалу
w = 0:Fs/FftL:Fs - Fs/FftL; % "вырезаем" реальные частоты >= 0

% Параметры передаточной функции скважины
    % коэффициент
    a = 0.4; % 1/с
    % глубина
    L = 1000; % м
    % модуль
    K = 4*10^(-8); % м^3/(Па*с)
    K_l = 2.2 * 10^9; % Па
    % модуль Юнга для материала трубы
    E = 2 * 10^11; % Па
    % постоянная времени компенсатора
    T = 0.01; % 1/с
    % Гидравлические сопротивления
     R = 1e6; % Па/(м^3/с)^2
      R_1 = R/2;
      R_L = R/2;
    % плотность бурового раствора
    ro = 2700; % кг/м^3
    K_pr = K_l / (1 + a*K_l/E);
    c = sqrt(K_pr / ro); % скорость звука в жидкости, м/с
    % диаметр трубы
    d = 114; % мм
    % площадь сечения трубы
    f = pi*d^2/4 *10^(-6); % м^2
    
    % Расчёт передаточной ф-ции
    v = 1/c*(w/2).^0.5 .* ((sqrt(w.^2+4*a^2)-w).^0.5 + 1i*(sqrt(w.^2+4*a^2)+w).^0.5);

    ro_g = c/f * (1 + 4*a^2 ./ w.^2).^0.25 .* exp(-1i*0.5*atan(2*a./w));

    W = (1i * w * K) ./ (1 + 1i*w.*T);
    

    % по давлению
    W_line = ( (1+2*R*W) .* cosh(v*L) + (W*ro .* ro_g + ...
    2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);

    W_line(1) = 1; % ручная коррекция NaN на нулевой частоте

    % по расходу
%     W_line = W .* ( (1+2*R*W) .* cosh(v*L) + (W*ro.*ro_g + ...
%     2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);
%     W_line = W_line * 10^8 * 0.8;

   

  


% --------------Пробуем сделать полосовой фильтр--------------------------
% order = 64; % Для большого шума 2048
% cut_off = (0.5) / (Fs/2); % cut_off fr-cy in pi units (0.5 Hz)
% cut_off_2 = (2.5) / (Fs/2);
% h = fir1(order, [cut_off cut_off_2]);

% Фурье фильтра
% fft_h = fft(h, FftL);
% fft_h_A = abs(fft_h) * 2 ./ FftL;
% fft_h_A(1) = fft_h_A(1) / 2;


% --------------------------исследуемый сигнал----------------------------
% Нужно выбрать тип сигнала
Signal = sin(2*pi*1.7*t);
% Signal = 5*sin(2*pi*1.7*t) + 5*sin(2*pi*10*t) + 2*sin(2*pi*0.7*t);

% ? how ? 5) Зададим фазовую модуляцию 0 1 0 на отрезке 30 сек
% Signal = [];
% Signal(1:length(t)/3) = 7*sin(2*pi*t(1:length(t)/3));
% Signal(length(t)/3+1:2*length(t)/3) = 7*sin(2*pi*t(length(t)/3+1:2*length(t)/3) - pi/4);
% Signal(2*length(t)/3 + 1:length(t)) = 7*sin(2*pi*t(2*length(t)/3 + 1:length(t)));

% 6) Зададим частотную модуляцию 0 1 0
% Signal = [];
% Signal(1:length(t)/3) = 7*sin(2*pi*t(1:length(t)/3));
% Signal(length(t)/3+1:2*length(t)/3) = 5*sin(2*pi*4*t(length(t)/3+1:2*length(t)/3));
% Signal(2*length(t)/3 + 1:length(t)) = 7*sin(2*pi*t(2*length(t)/3 + 1:length(t)));


% 1) простые прямоугольные импульсы
% T_s = 2;
% tau = 1;
% Signal = rem(t,T_s)<tau; % несимметричный импульс

% 2) половина в "1", половина в "0"
% Signal(1:length(t)/2) = rem(t(1:length(t)/2),T_s)<tau; % несимметричный импульс
% Signal(length(t)/2+1:length(t)) = zeros([1 length(t)/2]);

% 3) Зуб по середине
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
%находим вступление широкоимпульсной последовательности и конец
% [row, col] = find(t>=3+tau/2 & t<=23+tau/2);
% Signal(col(1):col(end)) = (rem(t(col(1):col(end)),T)<tau);
% Signal(1:col(1)-1) = 0;
% % Кусок коротких импульсов
% [row,col] = find(t>=23+tau/2 & t<=26+tau/2);
% Signal(col(1):col(end)) = (rem(t(col(1):col(end)),T/2)<tau/2);
% Signal(col(end)+1:length(t)) = 0;

% Образ сигнала идеальный (без шума)
original_fft = fft(Signal, FftL);

% Образ выходного сигнала после прохождения линии
Fft_S_out = W_line .* (original_fft);

% Входной по выходному сигналу с помощью передаточной ф-ции W_line (шума нет)
fft_restored = Fft_S_out ./ W_line;
restored_signal = ifft(fft_restored, FftL);

% Output signal in time domain 
S_out = ifft(Fft_S_out);
%----------------------------
% Создаём шум
err = rand(1,length(t));
err = err - mean(err);
err = err * real(max(S_out))*0.2 / std(err);
% Добавляем шум к выходному сигналу
S_out_noise = S_out + err;
% Преобразуем обратно в частотный домен
Fft_S_out_noise = fft(S_out_noise, FftL);


%-------------------------------filter apply------------------------------
% Fft_S_out_noise = Fft_S_out_noise .* fft_h;
%-------------------------------------------------------------------------


% Восстанавливаем входной сигнал по шумному выходному 
% с помощью передаточной функции
fft_restored_signal = Fft_S_out_noise ./ W_line;
restored_signal_noise = ifft(fft_restored_signal, FftL);
restored_signal_noise = real(restored_signal_noise); % Избавляемся от мнимых дребезгов


%-------------------------сглаживающий фильтр---------------------------
% cubic_test_filt = sgolayfilt(restored_signal_noise, 3, 251);

%-------------Дополнительный сглаживающий MA-фильтр----------------------
% counts_per_frame = 50;
% coef25MA = ones(1, counts_per_frame)/counts_per_frame;
% MA_sig = filter(coef25MA, 1, restored_signal_noise);



% АЧХ для образа исходного сигнала
original_fft_A = abs(original_fft);
original_fft_A = 2*original_fft_A ./ FftL;
original_fft_A(1) = original_fft_A(1)/2;

% АЧХ для образа Выходного сигнала
Fft_out_A=abs(Fft_S_out);% Амплитуды преобразования Фурье сигнала
Fft_out_A=2*Fft_out_A./FftL;% Нормировка спектра по амплитуде
Fft_out_A(1)=Fft_out_A(1)/2;% Нормировка постоянной составляющей в спектре

% Частоты
F = Fs*(0:FftL/2-1)/FftL; % shifting


%------расхождение сигналов с увеличением амплитуды погрешности----------
% sig_diff = abs( real(restored_signal_noise) - Signal);
% sig_diff_squares = sig_diff .^2;
% sq_err = sqrt(sum(sig_diff_squares));

%-----------------------filter АЧХ-----------------------------------------
% fft_h_A = abs(fft_h);
% fft_h_A = 2*fft_h_A ./ FftL;
% fft_h_A(1) = fft_h_A(1)/2;


%-----------------------------------АЧХ--------------------------------
% figure('Name', 'АЧХ для входного и выходного сигналов', 'NumberTitle', 'off');
% subplot(2,1,1);
% plot(F,original_fft_A(1:length(F)));
% title('АЧХ входного сигнала');
% xlabel('F, Hz');
% ylabel('A');
% 
% subplot(2,1,2);
% plot(F, Fft_out_A(1:length(F)));
% title('АЧХ выходного сигнала');
% xlabel('F, Hz');
% ylabel('A');

%-----------------------------------АЧХ фильтра------------------------
% figure('Name', 'АЧХ полосного фильтра')
% plot(F, fft_h_A(1:length(F)));
% title('АЧХ полосного фильтра');
% xlabel('F, Hz');
% ylabel('A');


% Просто смотрим, что подали на вход ПГИ
figure('Name', 'Идеальный входной сигнал', 'NumberTitle', 'off');
% check for original signal
S_orig = ifft(original_fft);
plot(t,S_orig);
title('Входной сигнал');
xlabel('t, s');
ylabel('A');


% ----------------getting signals in time domain---------------------
figure('Name', 'Входные сигналы, восстановленные по чистому и шумному выходному', 'NumberTitle', 'off');
subplot(2,1,1);
% Входной сигнал, восстановленный по ЧИСТОМУ выходному с помощью W_line
plot(t, real(restored_signal));
title('Входной сигнал из чистого выходного');
xlabel('t, s');
ylabel('A');
% ylim([-0.2 1.2]);

% Входной сигнал, восстановленный по выходному ШУМНОМУ с помощью W_line
subplot(2,1,2);
plot(t, (restored_signal_noise));
title('Входной сигнал из шумного выходного');
xlabel('t, s');
ylabel('A');

%-----------------------обзор выходных "страшненьких" сигналов-----------
figure('Name', 'Выходные сигналы: чистый и с шумом', 'NumberTitle', 'off');
% Выходной "страшненький" сигнал на поверхности
subplot(2,1,1);
plot(t, real(S_out));
title('Чистый выходной сигнал');
xlabel('t, s');
ylabel('A');

% Выходной, да ещё и с шумом
subplot(2,1,2);
plot(t, real(S_out_noise));
title('Шумный выходной сигнал');
xlabel('t, s');
ylabel('A');


% %-----------------------среднеквадратичная ошибка-----------------------
% figure('Name', 'Среднеквадратичная ошибка', 'NumberTitle', 'off');
% plot(A_err, sq_err);
% title('Чистый выходной сигнал');
% xlabel('t, s');
% ylabel('A');

%--------------------Применение сглаживающих фильтров----------------------
% figure('Name', 'Сглаживание восстановленного выходного сигнала', 'NumberTitle', 'off');
% plot(t, cubic_test_filt);
% title('Шумный выходной сигнал после дополнительного сглаживания');
% xlabel('t, s');
% ylabel('A');
% 
% %------------------Дополнительное сглаживание MA-фильтром-----------------
figure('Name', 'Сглаживание восстановленного выходного сигнала при помощи MA-фильтра', 'NumberTitle', 'off');
plot(t, MA_sig);
title('Шумный выходной сигнал после дополнительного сглаживания при помощи MA-фильтра');
xlabel('t, s');
ylabel('A');