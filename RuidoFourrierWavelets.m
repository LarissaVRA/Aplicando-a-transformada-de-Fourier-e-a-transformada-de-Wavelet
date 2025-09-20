clc; clear; close all;

% === 1. Ler o áudio ===
filename = 'topicos.wav';   
[x, Fs] = audioread(filename);  

% Verificar se o arquivo foi encontrado
if isempty(x)
    error('Arquivo "%s" não encontrado! Verifique o nome e a pasta.', filename);
end

% Se for estéreo, converter para mono
if size(x, 2) == 2
    x = mean(x, 2);
    disp('Áudio convertido de estéreo para mono');
end

% Normalizar sinal
x = x / max(abs(x));

% === 2. Adicionar ruído branco gaussiano ===
SNR_dB = 5;  % Relação sinal-ruído em dB
Px = mean(x.^2);
SNR_linear = 10^(SNR_dB/10);
Pn = Px / SNR_linear;
noise = sqrt(Pn) * randn(size(x));
x_noisy = x + noise;

% === 3. Redução de ruído via FFT ===
N = length(x_noisy);
X = fft(x_noisy);
f = (0:N-1)*(Fs/N);
fc = 2000;  % frequência de corte em Hz

H = (f <= fc | f >= (Fs - fc))';
if size(H, 2) > 1
    H = H';
end

X_filtered = X .* H;
x_denoised_fft = real(ifft(X_filtered));

% === 4. REDUÇÃO DE RUÍDO VIA WAVELET (SEM TOOLBOX) ===
% Método alternativo usando DWT manual simples

% Função simples para DWT manual
function [approx, detail] = simple_dwt(signal)
    % Filtros Haar wavelet simples
    h = [1/sqrt(2), 1/sqrt(2)];  % Filtro passa-baixa
    g = [1/sqrt(2), -1/sqrt(2)]; % Filtro passa-alta
    
    % Aplicar convolução e downsampling
    approx = conv(signal, h, 'same');
    approx = approx(1:2:end);
    
    detail = conv(signal, g, 'same');
    detail = detail(1:2:end);
end

% Aplicar DWT manual em 3 níveis
level = 3;
signals = cell(level + 1, 1);
signals{1} = x_noisy;

for i = 1:level
    [approx, detail] = simple_dwt(signals{i});
    signals{i+1} = approx;
end

% Aplicar thresholding nos coeficientes de detalhe
threshold = 0.1 * max(abs(x_noisy));  % Threshold adaptativo

for i = 2:level + 1
    % Aplicar soft thresholding
    coeffs = signals{i};
    signs = sign(coeffs);
    coeffs_thresh = signs .* max(abs(coeffs) - threshold, 0);
    signals{i} = coeffs_thresh;
end

% Reconstrução manual (simplificada)
x_denoised_wavelet = signals{level + 1};
for i = level:-1:1
    % Interpolar e reconstruir (aproximação simplificada)
    x_denoised_wavelet = interp1(1:length(x_denoised_wavelet), x_denoised_wavelet, ...
                                linspace(1, length(x_denoised_wavelet), length(signals{i})), 'spline');
end

% Ajustar tamanho
if length(x_denoised_wavelet) > N
    x_denoised_wavelet = x_denoised_wavelet(1:N);
elseif length(x_denoised_wavelet) < N
    x_denoised_wavelet(end+1:N) = 0;
end

% === 5. ALTERNATIVA: Filtro passa-banda como simulação wavelet ===
% Esta é uma alternativa mais robusta sem toolbox
disp('Usando método alternativo de filtragem multibanda...');

% Criar filtro passa-banda para simular efeito wavelet
f_low = 100;  % Frequência baixa
f_high = 4000; % Frequência alta

% Filtro passa-banda
H_bp = (f >= f_low & f <= f_high) | (f >= (Fs - f_high) & f <= (Fs - f_low));
H_bp = H_bp';

X_filtered_bp = X .* H_bp;
x_denoised_wavelet = real(ifft(X_filtered_bp));

% === 6. Plotar comparações ===
t = (0:N-1)/Fs;

figure('Position', [100, 100, 1200, 800]);

subplot(4,1,1);
plot(t, x);
title('Sinal Original');
xlabel('Tempo [s]'); ylabel('Amplitude');

subplot(4,1,2);
plot(t, x_noisy);
title(sprintf('Sinal com Ruído Branco (SNR = %d dB)', SNR_dB));
xlabel('Tempo [s]'); ylabel('Amplitude');

subplot(4,1,3);
plot(t, x_denoised_fft);
title('Redução de Ruído - Método FFT (Passa-Baixa)');
xlabel('Tempo [s]'); ylabel('Amplitude');

subplot(4,1,4);
plot(t, x_denoised_wavelet);
title('Redução de Ruído - Método Multibanda (Alternativa Wavelet)');
xlabel('Tempo [s]'); ylabel('Amplitude');

% === 7. Calcular métricas de qualidade ===
fprintf('=== MÉTRICAS DE QUALIDADE ===\n');

snr_noisy = 10*log10(mean(x.^2)/mean((x_noisy-x).^2));
snr_fft = 10*log10(mean(x.^2)/mean((x_denoised_fft-x).^2));
snr_wavelet = 10*log10(mean(x.^2)/mean((x_denoised_wavelet-x).^2));

fprintf('SNR Ruidoso: %.2f dB\n', snr_noisy);
fprintf('SNR após FFT: %.2f dB\n', snr_fft);
fprintf('SNR após Multibanda: %.2f dB\n', snr_wavelet);

% === 8. Reproduzir áudios ===
disp(' ');
disp('=== REPRODUÇÃO DOS ÁUDIOS ===');

disp('Reproduzindo: Original');
sound(x, Fs);
pause(N/Fs + 1);

disp('Reproduzindo: Com Ruído');
sound(x_noisy, Fs);
pause(N/Fs + 1);

disp('Reproduzindo: Reduzido por FFT');
sound(x_denoised_fft, Fs);
pause(N/Fs + 1);

disp('Reproduzindo: Reduzido por Método Multibanda');
sound(x_denoised_wavelet, Fs);