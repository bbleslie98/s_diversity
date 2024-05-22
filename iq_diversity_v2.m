clearvars;
close all;

path = "..\IQs\Out\";
file_prefix = "09_04_4Mfs_200kbps_";
% ix_s = 3.6e6;
% ix_e = 6.7e6;
ix_s = 2e6;
ix_e = 9e6;

% file_prefix = "09_04_800ksps_200kbps_snr3db_";
% ix_s = 0.7e6;
% ix_e = 1.4e6;

%% Import data

% Signal details
fs = 4e6;
DR = 200e3;
spb = fs / DR;


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

%% Gen syncword IQ

s_spb = 2;

fgauss = gaussdesign(0.5, 2, s_spb);

syncstr     = '11100011000111001001110110101110';
% syncstr     = '1010101010101010101010101010101011100011000111001001110110101110';
slen = length(syncstr);

sync = zeros(1, slen);

for i = 1:slen 
    sync(i) = str2num(syncstr(i));
end
sync = sync*2-1;

sm = zeros(1,slen*s_spb);

int_sync = 0;
for k = 1:slen*s_spb

    int_sync = int_sync + sync(ceil(k/s_spb));
    sm(k) = int_sync;
end

fc = 500e3;

t = (0:slen*s_spb-1)/fs;
sync_bb = exp(1j* 2*pi*DR/4./fs.*sm);

% sync_bb = filter(fgauss, 1, sync_bb);

sync_c = sync_bb .* exp(1j* 2*pi*fc*t);

%% Process

% corr = xcorr(s_ch1, s_ch2);
% [mv, mi] = max(abs(corr));
% rot = angle(corr(mi));
s_ch2 = s_ch2 * exp( 1j*pi/180*80);
% % s_ch2 = awgn(s_ch2, -5, 30);
% 
% corr = xcorr(s_ch1, s_ch2);
% [mv, mi] = max(abs(corr));


freq_x = ( 1:length(s_ch1) ) * fs / len; 

% b = fir1(60, 700e3*2/fs);
% 
% s_ch1 = filter(b, 1, s_ch1);
% s_ch2 = filter(b, 1, s_ch2);

spect_1 = fft((s_ch1));
spect_2 = fft((s_ch2));

% figure;
% plot(freq_x, abs(spect_1));
% figure;
% plot(freq_x, abs(spect_2));


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

[x1, xi] = xcorr(bb2, sync_bb);
cmp = 3*movmean(abs(x1), 10000);

absval = abs(x1);

figure;
plot(xi, absval);
hold on;
plot(xi, cmp);

% compare levels
ov_cmp_ids =  (cmp > 200) & (cmp < absval);
% highval = (absval - cmp).*packet_ids;

figure;
stem(ov_cmp_ids);


% PACKET STREAM FILTER psf
% 2 sample per bit: 25 packet average: 4524,528 samples between packets

psf_pck_spacing = 4524;
psf_pck_peakwin = 31;
psf_peaknum     = 10;
psf_resol_win   = 6;
psf_overlap     = 2;
psf_ext_factor  = 0.1;  % Extension factor extends the end of the resolution interval by the given fraction
% Overlap looks for maximums from the previous interval
% Total resolution win is resol win + overlap

psf_pck_blankiv = psf_pck_spacing - psf_pck_peakwin;

filter = 0;

for i = 1:psf_peaknum
    filter = [filter, triang(psf_pck_peakwin)', zeros(1, psf_pck_blankiv) ];

end

figure;
[filtcorr, fitlcorr_id] = xcorr(ov_cmp_ids, filter);
filtcorr = abs(filtcorr);
plot(filtcorr);


valid_filt_id = abs(filtcorr) > 15;
nz_valid_filt_id = find(valid_filt_id); % non zero elements of the id

pk_sync_ids = [];
for i = nz_valid_filt_id(1) : psf_pck_spacing * psf_resol_win : nz_valid_filt_id(end)
    [mv, mi] = max( filtcorr( (i - ceil(psf_pck_spacing * (psf_overlap + psf_ext_factor))) : ( i + ceil(psf_pck_spacing * (psf_resol_win + psf_ext_factor))) ) );
    pk_sync_ids = [pk_sync_ids, i + mi];
end

pk_sync_padding = 4;
% First ID is when the filter value starts to be valid
pk_sync_ids = [ pk_sync_ids(1) - ceil( (pk_sync_ids(1) - nz_valid_filt_id(1))/ psf_pck_spacing + pk_sync_padding) * psf_pck_spacing, pk_sync_ids];
% Appending the last id, when filter value is still valid
pk_sync_ids = [ pk_sync_ids, pk_sync_ids(end) + ceil((nz_valid_filt_id(end) - pk_sync_ids(end))/psf_pck_spacing + pk_sync_padding) * psf_pck_spacing];

pk_ids = [];
for i = 1:length(pk_sync_ids)-1
   
    j = 0;
    while pk_sync_ids(i) + j * psf_pck_spacing < pk_sync_ids(i+1)

        pk_ids = [pk_ids, pk_sync_ids(i) + j * psf_pck_spacing];
        j = j+1;
    end
end



q = zeros(1, pk_sync_ids(end));
q(pk_ids) = 1;
% hold on;
% plot(100*q, "LineStyle","none","Marker","+");
figure;
stem(q)

% PKT indexes transformed back to the bb indexes
pkt_id_norm = xi(fitlcorr_id(pk_ids));

% Remove packets which are not in the comparation range
p = find(ov_cmp_ids);   % Indexes which are over the cmp value
% xi(p(1));
% xi(p(end));

pkt_id_norm = pkt_id_norm(pkt_id_norm > xi(p(1)) & pkt_id_norm < xi(p(end)));

figure;
plot(real(bb2));
hold on;
plot(pkt_id_norm, zeros(1, length(pkt_id_norm)), "LineStyle","none", "Marker","+","MarkerSize",8);

% finding the packet power levels
% s_ch1(10*(pkt_id_norm(1)+(-420:(4160-400))))

pck_count = length(pkt_id_norm);
pwr1s = zeros(pck_count,1);
pwr2s = pwr1s;
phases = pwr1s;

tcnco = ( 0 : (len-1) )/fs;
if1_fs = s_ch1 .* exp(1j*( fs - f_if ) .* tcnco');
if2_fs = s_ch2 .* exp(1j*( fs - f_if ) .* tcnco');

for i = 1:pck_count
    pkt1 = s_ch1(10*(pkt_id_norm(i)+(-420:(4160-401))));
    pkt2 = s_ch2(10*(pkt_id_norm(i)+(-420:(4160-401))));
    
    pwr1s(i) = pkt1'*pkt1 / length(pkt1);
    pwr2s(i) = pkt2'*pkt2 / length(pkt2);

    corr = xcorr(pkt1, pkt2);
    [mv, mi] = max(abs(corr));

    phases(i) = 180/pi * angle( corr(mi) );

    % FOR NOW THE DIFFERENCES BETWEEN 500KHZ IF AND FS IF ANGLES IS
    % NUMERICAL ONLY
%     pkt1 = if1_fs(10*(pkt_id_norm(i)+(-420:(4160-401))));
%     pkt2 = if2_fs(10*(pkt_id_norm(i)+(-420:(4160-401))));
% 
%     corr = xcorr(pkt1, pkt2);
%     [mv, mi] = max(abs(corr));
% 
%     phases(i) = phases(i) - 180/pi * angle( corr(mi) );
end

figure
stem(phases)


% plot(3*movmean(abs(x1), 5000) > abs(x1));

% figure;
% plot(demod1(1:2:end));
% figure;
% plot(demod1(2:2:end)),
% 
% [scorr, scorr_lags] = xcorr(demod1(1:2:end-1), sync);
% [scorr(2, :), scorr_lags(2, :)] = xcorr(demod1(2:2:end), sync);
% [scorr(3, :), scorr_lags(3, :)] = xcorr(demod2(1:2:end-1), sync);
% [scorr(4, :), scorr_lags(4, :)] = xcorr(demod2(2:2:end), sync);
% 
% tiledlayout(2,2);
% 
% for i = 1:4
%     nexttile
%     plot(scorr_lags(i,:), abs(scorr(i,:)));
% end


