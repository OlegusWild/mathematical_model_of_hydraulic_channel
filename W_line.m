% ??????? ??? ???????????? ???????????? ?-??? ??????????????? ?????? ?????
% ?? ???? ????????: 
% 1) ????????? ????????: ??????? L[?], ????????? ???????? ???????? ro[??/?^3], ???????
% ???????????? d[??]
% 2) ????????? ?????-??????????????: ??????? ????????????? Fs[], ?????????
% ???????? duration[beg,? end,?]
% 3) ?????: type - 'P' ??? 'Q' (?? ???????? ??? ?? ???????)

% ?????? ?????????????: 
% [W_line, w, t] = W_line(1000, 2700, 114, 1000, [0 10], 'P');

function [y, w, t, FftL]= W_line(L, ro, d, Fs, duration, type)

    Ts = 1/Fs; % ??? ?? ???????
    FftL = 2^nextpow2(length(duration(1):Ts:duration(2))); % ???-?? ????? ??????????????
    t = linspace(duration(1),duration(2), FftL); % ? ???????????? ? ???-??? ????? ?????? ????????? ?????
    w = 0:Fs/FftL:Fs - Fs/FftL; % "????????" ???????? ??????? >= 0

    % ?????? ????????? ???????????? ??????? ????????, ??????? ??????? ????????
    % ??? ?????????????? ?? ??????? ??????????
        % ??????????? ?????????
        a = 0.4; % 1/?
		% ?????? ????????? ?????? ????????
        K_l = 2.2 * 10^9; % ??
        % ?????? ???? ??? ????????? ????? (?????)
        E = 2 * 10^11; % ??
		% ???????? ????????????
        K = 4*10^(-8); % ?^3/(??*?)
        % ?????????? ??????? ????????????
        T = 0.01; % 1/?
        % ?????????????? ?????????????
         R = 1e6; % ??/(?^3/?)^2
          R_1 = R/2;
          R_L = R/2;
		
        K_pr = K_l / (1 + a*K_l/E);
        c = sqrt(K_pr / ro); % ???????? ????? ? ????????, ?/?
        % ??????? ??????? ?????
        f = pi*d^2/4 *10^(-6); % ?^2

     % ?????? ???????????? ?-???
        v = 1/c*(w/2).^0.5 .* ((sqrt(w.^2+4*a^2)-w).^0.5 + 1i*(sqrt(w.^2+4*a^2)+w).^0.5);

        ro_g = c/f * (1 + 4*a^2 ./ w.^2).^0.25 .* exp(-1i*0.5*atan(2*a./w));

        W = (1i * w * K) ./ (1 + 1i*w.*T);

    if type == 'P'
        % ?? ????????
        W_line = ( (1+2*R*W) .* cosh(v*L) + (W*ro .* ro_g + ...
        2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);

        W_line(1) = 1; % ?????? ????????? NaN ?? ??????? ???????
    end

    if type == 'Q'
        % ?? ???????
        W_line = W .* ( (1+2*R*W) .* cosh(v*L) + (W*ro.*ro_g + ...
        2*R_L ./ (ro*ro_g) .* (1+2*R_1*W) ) .* sinh(v*L)).^(-1);
        % W_line = W_line * 10^8 * 0.8; ? ??? ?????????????? ? ????????????
    end
    
    y = W_line;
