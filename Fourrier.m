clear
clc

A=imread('teste.jpg');
B=double(rgb2gray(A)); % Convert RBG->gray, 256 bit->double.
Bt=fft2(B); % B is grayscale image from above
Btsort = sort(abs(Bt(:))); % Sort by magnitude
% Zero out all small coefficients and inverse transform
for keep=[.1 .05 .01 .002];
thresh = Btsort(floor((1-keep)*length(Btsort)));
ind = abs(Bt)>thresh; % Find small indices
Atlow = Bt.*ind; % Threshold small indices
Alow=uint8(ifft2(Atlow)); % Compressed image
figure, imshow(Alow) % Plot Reconstruction
end

