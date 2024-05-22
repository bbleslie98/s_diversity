clearvars;
close all;

path = "..\IQs\Out\";
file_prefix = "09_04_4Mfs_200kbps_";
ix_s = 3.6e6;
ix_e = 6.7e6;

% file_prefix = "09_04_800ksps_200kbps_snr3db_";
% ix_s = 0.7e6;
% ix_e = 1.4e6;

%% Import data

file_ch1 = path + file_prefix + "ch0.cf32";
file_ch2 = path + file_prefix + "ch1.cf32";

fch1 = fopen(file_ch1);
fch2 = fopen(file_ch2);

data_ch1 = fread(fch1, 'float32');
data_ch2 = fread(fch2, 'float32');

s_ch1 = data_ch1(1:2:end) + data_ch1(2:2:end) * 1j;
s_ch2 = data_ch2(1:2:end) + data_ch2(2:2:end) * 1j;

%% Cutout

s_ch1 = s_ch1(ix_s:ix_e);
s_ch2 = s_ch2(ix_s:ix_e);

len = length(s_ch1);

%% Process

fs = 4e6;
DR = 200e3;
spb = fs / DR;

corr = xcorr(s_ch1, s_ch2);

[mv, mi] = max(abs(corr));

rot = angle(corr(mi));

s_ch2 = s_ch2 * exp( 1j*rot);

corr = xcorr(s_ch1, s_ch2);
[mv, mi] = max(abs(corr));


freq_x = ( 1:length(s_ch1) ) * fs / len; 

% b = fir1(60, 700e3*2/fs);
% 
% s_ch1 = filter(b, 1, s_ch1);
% s_ch2 = filter(b, 1, s_ch2);

spect_1 = fft((s_ch1));
spect_2 = fft((s_ch2));

figure;
plot(freq_x, abs(spect_1));
figure;
plot(freq_x, abs(spect_2));


%% Dec

if1 = s_ch1;
if2 = s_ch2;

t_spb = 2;
dec = spb / t_spb;
declen = floor(len/dec);

ifd1 = zeros(1, declen);
ifd2 = ifd1;

for i = 1:declen
    
    if i*dec > len
        break;
    end
    ifd1(i) = sum(if1( (i-1)*dec + 1 : i*dec ))/dec;
    ifd2(i) = sum(if2( (i-1)*dec + 1 : i*dec ))/dec;
end

%% BB conv

f_if = 500e3;

tcnco = (0:declen-1)/fs/dec;

bb1 = ifd1 .* exp(-1j*2*pi*f_if.*tcnco);
bb2 = ifd2 .* exp(-1j*2*pi*f_if.*tcnco);


%% Demod

syncstr = '11100011000111001001110110101110';
sync = zeros(1, length(syncstr));
for i = 1:length(syncstr) 
    sync(i) = str2num(syncstr(i));
end
sync = sync*2-1;

demod1 = sign(imag(bb1(1:end-1).*conj(bb1(2:end))));
demod2 = sign(imag(bb2(1:end-1).*conj(bb2(2:end))));

figure;
plot(demod1(1:2:end));
figure;
plot(demod1(2:2:end)),

[scorr, scorr_lags] = xcorr(demod1(1:2:end-1), sync);
[scorr(2, :), scorr_lags(2, :)] = xcorr(demod1(2:2:end), sync);
[scorr(3, :), scorr_lags(3, :)] = xcorr(demod2(1:2:end-1), sync);
[scorr(4, :), scorr_lags(4, :)] = xcorr(demod2(2:2:end), sync);

tiledlayout(2,2);

for i = 1:4
    nexttile
    plot(scorr_lags(i,:), abs(scorr(i,:)));
end
