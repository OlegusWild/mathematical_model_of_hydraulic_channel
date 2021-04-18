% Типы исследуемых сигналов
% На вход подаётся название сигнала (см. описание ниже в ф-ции)
% Ф-ция возвращает сигнал на заданной временной шкале

function Signal = signal_type(t, sig_name)
if sig_name == "sine"
    % 1) Просто синус
    Signal = sin(2*pi*1.7*t);
end

if sig_name == "harmonical"
    % 2) Сигнал с несколькими гармониками
    Signal = 5*sin(2*pi*1.7*t) + 5*sin(2*pi*10*t) + 2*sin(2*pi*0.7*t);
end

if sig_name == "rectangles"
    % 3.1) Простые прямоугольные импульсы
    T_s = 2;
    tau = 1;
    Signal = rem(t,T_s)<tau; % несимметричный импульс
end

if sig_name == "one_zero"
    % 3.2) половина в "1", половина в "0"
    Signal(1:length(t)/2) = rem(t(1:length(t)/2),T_s)<tau; % несимметричный импульс
    Signal(length(t)/2+1:length(t)) = zeros([1 length(t)/2]);
end

if sig_name == "tooth"
    % 3.3) Зуб по середине
    tau = 0.2; % длительность зуба
    pos = 5; % место на временной шкале, где расположен центр импульса
    [row, col] = find(t<pos+tau/2 & t>pos-tau/2);
    Signal(col(1):col(end)) = 1; % tooth
    Signal(1:col(1)-1) = 0;
    Signal(col(end)+1:length(t)) = 0;
end

if sig_name == "synch"
    % 3.4) SyncSequence_1
    tau = 2;
    T = 2*tau;
    % duration 30 s
    %находим вступление широкоимпульсной последовательности и конец
    [row, col] = find(t>=3+tau/2 & t<=23+tau/2);
    Signal(col(1):col(end)) = (rem(t(col(1):col(end)),T)<tau);
    Signal(1:col(1)-1) = 0;
    % % Кусок коротких импульсов
    [row,col] = find(t>=23+tau/2 & t<=26+tau/2);
    Signal(col(1):col(end)) = (rem(t(col(1):col(end)),T/2)<tau/2);
    Signal(col(end)+1:length(t)) = 0;
end





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




