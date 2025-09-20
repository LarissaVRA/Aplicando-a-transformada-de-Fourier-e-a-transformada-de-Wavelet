clear
clc

% Carregar e preparar a imagem
A = imread('teste.jpg');
B = double(rgb2gray(A)); % Convert RBG->gray, 256 bit->double

% Parâmetros
n_levels = 2; % Reduzindo para 2 níveis para melhor estabilidade
keep_values = [1.0 0.1 0.01 0.001]; % Valores mais altos para testar

% Aplicar compressão wavelet para diferentes taxas
for keep_idx = 1:length(keep_values)
    keep = keep_values(keep_idx);
    
    % Aplicar transformada wavelet simplificada (Haar)
    [wavelet_coeffs, sizes] = wavelet_decomposition(B, n_levels);
    
    % Achatar e ordenar coeficientes por magnitude
    all_coeffs = wavelet_coeffs(:);
    Csort = sort(abs(all_coeffs));
    
    % Calcular threshold
    thresh = Csort(max(1, floor((1-keep)*length(Csort))));
    
    % Aplicar threshold - zerar coeficientes pequenos
    wavelet_coeffs_compressed = wavelet_coeffs .* (abs(wavelet_coeffs) >= thresh);
    
    % Reconstruir a imagem
    Alow = wavelet_reconstruction(wavelet_coeffs_compressed, sizes, n_levels);
    
    % Limitar valores ao intervalo válido e converter
    Alow = max(0, min(255, Alow));
    Alow_uint8 = uint8(Alow);
    
    % Mostrar a imagem reconstruída
    figure
    imshow(Alow_uint8)
    title(['Compressão com ', num2str(keep*100), '% dos coeficientes wavelet'])
    
    % Calcular e mostrar métrica de qualidade
    mse = mean((B(:) - Alow(:)).^2);
    psnr = 10 * log10(255^2 / mse);
    fprintf('Para %.1f%%: MSE = %.2f, PSNR = %.2f dB\n', keep*100, mse, psnr);
end

% Função de decomposição wavelet melhorada
function [coeffs, sizes] = wavelet_decomposition(img, n_levels)
    coeffs = img;
    sizes = zeros(n_levels, 2);
    
    for level = 1:n_levels
        [rows, cols] = size(coeffs);
        sizes(level, :) = [rows, cols];
        
        % Aplicar wavelet Haar por linhas
        temp = zeros(rows, cols);
        for i = 1:rows
            temp(i, :) = haar_1d_forward(coeffs(i, :));
        end
        
        % Aplicar wavelet Haar por colunas
        for j = 1:cols
            temp(:, j) = haar_1d_forward(temp(:, j)')';
        end
        
        % Reorganizar os coeficientes
        half_r = ceil(rows/2);
        half_c = ceil(cols/2);
        
        LL = temp(1:half_r, 1:half_c);
        LH = temp(1:half_r, half_c+1:cols);
        HL = temp(half_r+1:rows, 1:half_c);
        HH = temp(half_r+1:rows, half_c+1:cols);
        
        % Combinar em uma única matriz
        coeffs = [LL, LH; HL, HH];
    end
end

% Função de reconstrução wavelet melhorada
function img = wavelet_reconstruction(coeffs, sizes, n_levels)
    img = coeffs;
    
    for level = n_levels:-1:1
        [rows, cols] = size(img);
        orig_rows = sizes(level, 1);
        orig_cols = sizes(level, 2);
        
        half_r = ceil(orig_rows/2);
        half_c = ceil(orig_cols/2);
        
        % Extrair sub-bandas
        LL = img(1:half_r, 1:half_c);
        LH = img(1:half_r, half_c+1:2*half_c);
        HL = img(half_r+1:2*half_r, 1:half_c);
        HH = img(half_r+1:2*half_r, half_c+1:2*half_c);
        
        % Reconstruir matriz completa
        reconstructed = zeros(orig_rows, orig_cols);
        reconstructed(1:half_r, 1:half_c) = LL;
        reconstructed(1:half_r, half_c+1:orig_cols) = LH(1:half_r, 1:orig_cols-half_c);
        reconstructed(half_r+1:orig_rows, 1:half_c) = HL(1:orig_rows-half_r, 1:half_c);
        reconstructed(half_r+1:orig_rows, half_c+1:orig_cols) = HH(1:orig_rows-half_r, 1:orig_cols-half_c);
        
        % Aplicar inversa por colunas
        for j = 1:orig_cols
            reconstructed(:, j) = haar_1d_inverse(reconstructed(:, j));
        end
        
        % Aplicar inversa por linhas
        for i = 1:orig_rows
            reconstructed(i, :) = haar_1d_inverse(reconstructed(i, :));
        end
        
        img = reconstructed;
    end
end

% Transformada Haar 1D direta
function result = haar_1d_forward(signal)
    n = length(signal);
    result = zeros(1, n);
    
    for i = 1:floor(n/2)
        result(i) = (signal(2*i-1) + signal(2*i)) / sqrt(2);      % Média (aproximação)
        if 2*i <= n
            result(i + floor(n/2)) = (signal(2*i-1) - signal(2*i)) / sqrt(2); % Diferença (detalhe)
        end
    end
    
    % Se n é ímpar, manter o último elemento
    if mod(n, 2) == 1
        result(floor(n/2) + 1) = signal(n) / sqrt(2);
    end
end

% Transformada Haar 1D inversa
function result = haar_1d_inverse(signal)
    n = length(signal);
    result = zeros(1, n);
    half = ceil(n/2);
    
    for i = 1:floor(n/2)
        result(2*i-1) = (signal(i) + signal(i + half)) / sqrt(2);
        result(2*i) = (signal(i) - signal(i + half)) / sqrt(2);
    end
    
    % Se n é ímpar, processar o último elemento
    if mod(n, 2) == 1
        result(n) = signal(half) * sqrt(2);
    end
end