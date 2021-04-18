% Сглаживающий фильтр "скользящее среднее"
% Принимает сигнал, который необходимо сгладить, а также размер окна

function smooth_sig = MA_filter(signal, counts_per_frame)
    coef25MA = ones(1, counts_per_frame)/counts_per_frame;
    smooth_sig = filter(coef25MA, 1, signal);