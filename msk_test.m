clearvars;
close all;

% f_if / fs * 360 => phase error step resolution

fs = 4e6;
DR = 200e3;

fc = 2267e3;

spb = fs/DR;

% preamble    = "10101010101010101010101010101010";
syncstr     = "1010101010101010101010101010101011100011000111001001110110101110";

packet      = "000011101010101111100110010000110010100110101000011001100000000101000101100100110110110110001101111010101001110001000100010111010111100011010010001010001010000101000010011110000111011011010101100011110101010101111011001000111011111010100101111001001001011101101110100000010101111000010001001111010000001001010000000000100100101011111110110011011111011110111110000101101000100001111011011011001000001001000100111010100010011011000111001000101110101011011111010000001011011100001101100000110000111001111110100001010010110110110010010110100000110010101111110101011100111110010001001000111110111010011111000001100110011111100010001010111000011111001011011111000010110001111110110101101110000110110011001001100010010000010000100001001000011000100100110110011001100101100001011000111001010001101011010001101000100010010010010110001010111011101101101001100110001001011011111011001000100011100111010000100111010100010010110111000100110001100000001001001101111101010001001101101011101000101010111100000100110010101000100001011001110100010100010100110010111001100001111000011110010110110010010111111111100110001010011001001000011101110101010001000101010101111111111011111000011010100011100011011011101010111010000100100101010001111111101101000000011000010110111000010010000101011100001111010111110000110010010111001101011010000000011110010011001011110100100001011000000111000110001000000111101011111001001100010010001110100100111000010111111010001001001101111010100100010110101001000011100111001100110111111110100010100011010011010011000001100110011001010101010110110100011011110110000010110111111101011100100011000010101111110011010001111000101000000100000001010011101100101010101100001000000001000011101100010110111100100100010100001000100010101101111001001001100001110011011100101011100010111110001111110011100110000001000001001111010101111111001101010100100000100010011101000110010101011100000010011110111000011001110001000110101110000111101111101010011001001100100110100011100101100011000001110100111100100011111000011011";

codestr = char(syncstr + packet);
syncstr = char(syncstr);

% codestr = '10111';

clen = length(codestr);
slen = length(syncstr);

code = zeros(1, clen);
sync = zeros(1, slen);

for i = 1:clen
    code(i) = str2num(codestr(i));
end
code = code*2-1;

for i = 1:slen 
    sync(i) = str2num(syncstr(i));
end
sync = sync*2-1;

sm = zeros(1,clen*spb);

int_code = 0;
for k = 1:clen*spb

    int_code = int_code + code(ceil(k/spb));
    sm(k) = int_code;
end

t = (0:clen*spb-1)/fs;
s_bb = exp(1j* 2*pi*DR/4./fs.*sm);
s_c = s_bb .* exp(1j* 2*pi*fc*t);

% plot(t, real(s_c*exp(1j*0*pi/180)));

%% Channel

% r_c = awgn(s_c, 30, 0);
snr1 = -5;
snr2 = -5;

N = 100;
angles = zeros(1,N);
% demodres = zeros(6,N);

for it = 1:N

    r_c1 = s_c  * 10^(snr1/20);
    r_c2 = s_c*exp(1j*0*pi/180) * 10^(snr2/20);
    
    % Bad, cuz the noise is not constant but the signal is
%     r_c1 = awgn(r_c1, snr1, 0); 
%     r_c2 = awgn(r_c2, snr2, 0);
    
    % 0 dBm power noise
    r_c1 = r_c1 + wgn(1, clen*spb, 0, 'complex');
    r_c2 = r_c2 + wgn(1, clen*spb, 0, 'complex');
    
    %% Demod
    
    f_if = 4000e3;
    
    f_lo = fc - f_if;
    
    r_if1 = r_c1 .* exp(-1j*2*pi*f_lo.*t);
    r_if2 = r_c2 .* exp(-1j*2*pi*f_lo.*t);
    
% correlation without decimation
    corr = xcorr(r_if1, r_if2);
    
    [mv, mi] = max(abs(corr));
%     disp( angle(corr(mi)) * 180/pi);
    angles(it) = angle(corr(mi)) * 180/pi;

end
plot(angles, LineStyle="none", Marker="x");
grid on;
figure;
histogram(angles, 512);
title("Avg: " + num2str(sum(angles)/N));
% end of correlation without decimation

% Dec
% Target samples per bit
t_spb = 2;
dec = floor(spb/t_spb);

r_ifd1 = zeros(1, floor(clen*t_spb));
r_ifd2 = r_ifd1;

for i = 1:clen*t_spb
    r_ifd1(i) = sum(r_if1( (i-1)*dec + 1 : i*dec ))/dec;
    r_ifd2(i) = sum(r_if2( (i-1)*dec + 1 : i*dec ))/dec;
end


% corerlation with decimation
%     corr = xcorr(r_ifd1, r_ifd2);
%     
%     [mv, mi] = max(abs(corr));
% %     disp( angle(corr(mi)) * 180/pi);
%     angles(it) = angle(corr(mi)) * 180/pi;
% 
% end
% stem(angles);
% figure;
% histogram(angles, 256);
% 
% return;
% end of correlation with decimation

tdec = (0:clen*t_spb-1)/fs*dec;

r_bb1 = r_ifd1 .* exp(-1j*2*pi*f_if.*tdec);
r_bb2 = r_ifd2 .* exp(-1j*2*pi*f_if.*tdec);

r_ifdp = r_ifd1 + r_ifd2 * exp(-1j*0*pi/180) / (10^((snr1-snr2)/20));

% /20 : 524.82
% /18 : 524.599
% /22 : 525.786

r_bbp = r_ifdp .* exp(-1j*2*pi*f_if.*tdec);


% plot(abs(fft(s_bb)));
% figure;
% plot(abs(fft(r_bb)));
% 
% figure;
% plot(real(r_bb))
% figure;
% plot(real(s_bb))

demod1 = sign(imag(r_bb1(1:end-1).*conj(r_bb1(2:end))));
demod2 = sign(imag(r_bb2(1:end-1).*conj(r_bb2(2:end))));

demodp = sign(imag(r_bbp(1:end-1).*conj(r_bbp(2:end))));

% iterated evaluation
% demodres(1,it) =  sum(abs( demodp(1:2:end-1) + code(1:end-1) ));
% demodres(2,it) =  sum(abs( demodp(2:2:end) + code(1:end-1) ));
% demodres(3,it) =  sum(abs( demod1(1:2:end-1) + code(1:end-1) ));
% demodres(4,it) =  sum(abs( demod1(2:2:end) + code(1:end-1) ));
% demodres(5,it) =  sum(abs( demod2(1:2:end-1) + code(1:end-1) ));
% demodres(6,it) =  sum(abs( demod2(2:2:end) + code(1:end-1) ));
% 
% end
% 
% figure;
% 
% titlelist = ["div a ", "div b ", "pol1 a ", "pol1 b ", "pol2 a ", "pol2 b "];
% 
% tiledlayout(3,2)
% 
% for i = 1:6
%     nexttile
%     plot(demodres(i,:));
%     title(titlelist(i) + num2str(sum(demodres(i,:))/N))
% end
% end of iterated evaluation

% subplot(2, 3, 1);
% plot(demodres(1,:));
% title("div a " + num2str(sum(demodres(1,:))/N));
% 
% subplot(2, 3, 4);
% plot(demodres(2,:));
% title("div b " + num2str(sum(demodres(2,:))/N));
% 
% subplot(2, 3, 2);
% plot(demodres(3,:));
% title("pol1 a " + num2str(sum(demodres(3,:))/N));
% 
% subplot(2, 3, 5);
% plot(demodres(4,:));
% title("pol1 b " + num2str(sum(demodres(4,:))/N));
% 
% subplot(2, 3, 3);
% plot(demodres(5,:));
% title("pol2 a " + num2str(sum(demodres(5,:))/N));
% 
% subplot(2, 3, 6);
% plot(demodres(6,:));
% title("pol2 b " + num2str(sum(demodres(6,:))/N));

% figure;
% subplot(2, 3, 1);
% cerrv = demodp(1:2:end-1) + code(1:end-1);
% stem(cerrv);
% title("pol a " + num2str(sum(abs(cerrv))));
% 
% subplot(2, 3, 4);
% cerrv = demodp(2:2:end) + code(1:end-1);
% stem(cerrv);
% title("pol b " + num2str(sum(abs(cerrv))));
% 
% subplot(2, 3, 2);
% cerrv = demod1(1:2:end-1) + code(1:end-1);
% stem(cerrv);
% title("dem1 a " + num2str(sum(abs(cerrv))));
% 
% subplot(2, 3, 5);
% cerrv = demod1(2:2:end) + code(1:end-1);
% stem(cerrv);
% title("dem1 b " + num2str(sum(abs(cerrv))));
% 
% subplot(2, 3, 3);
% cerrv = demod2(1:2:end-1) + code(1:end-1);
% stem(cerrv);
% title("dem2 a " + num2str(sum(abs(cerrv))));
% 
% subplot(2, 3, 6);
% cerrv = demod2(2:2:end) + code(1:end-1);
% stem(cerrv);
% title("dem2 b " + num2str(sum(abs(cerrv))));
