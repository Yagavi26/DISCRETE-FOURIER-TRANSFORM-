# EXP 1 A : COMPUTATION OF DFT USING DIRECT AND FFT

# AIM: 
To Obtain DFT and FFT of a given sequence in SCILAB.

# To Obtain DFT and FFT of a given sequence in SCILAB. 

# APPARATUS REQUIRED: 
PC installed with SCILAB. 

# PROGRAM: 
// DISCRETE FOURIER TRANSFORM 
clc;
clear;

xn = input("Enter the sequence in square brackets [ ] : ");

n1 = 0:1:length(xn)-1;


subplot(3,1,1);
plot2d3(n1, xn);
xlabel('Time n');
ylabel('Amplitude xn');
title('Input Sequence');

j = %i;
N = length(xn);
Xk = zeros(1, N);

for k = 0:N-1
    for n = 0:N-1
        Xk(k+1) = Xk(k+1) + xn(n+1) * exp((-j*2*%pi*k*n)/N);
    end
end

disp("DFT Values (Xk):");
for k = 1:N
    realPart = real(Xk(k));
    imagPart = imag(Xk(k));
    
    if abs(imagPart) < 1e-10 then
        printf("%d", round(realPart));
    elseif imagPart > 0 then
        printf("%d+%dj", round(realPart), round(imagPart));
    else
        printf("%d%dj", round(realPart), round(imagPart));
    end
    
    if k < N then
        printf(" , ");
    else
        printf("\n");
    end
end

K1 = 0:1:length(Xk)-1;

magnitude = abs(Xk);
subplot(3,1,2);
plot2d3(K1, magnitude);
xlabel('frequency (Hz)');
ylabel('magnitude (gain)');
title('Magnitude Spectrum');

angle = atan(imag(Xk), real(Xk));
subplot(3,1,3);
plot2d3(K1, angle);
xlabel('frequency (Hz)');
ylabel('Phase');
title('Phase Spectrum');
DFT USING FFT

clc;
clear;

xn = input("Enter the sequence in square brackets [ ] : ");

n1 = 0:1:length(xn)-1;


subplot(3,1,1);
plot2d3(n1, xn);
xlabel('Time n');
ylabel('Amplitude xn');
title('Input Sequence');

j = %i;
N = length(xn);
Xk = zeros(1, N);

for k = 0:N-1
    for n = 0:N-1
        Xk(k+1) = Xk(k+1) + xn(n+1) * exp((-j*2*%pi*k*n)/N);
    end
end

disp("DFT Values (Xk):");
for k = 1:N
    realPart = real(Xk(k));
    imagPart = imag(Xk(k));
    
    if abs(imagPart) < 1e-10 then
        printf("%d", round(realPart));
    elseif imagPart > 0 then
        printf("%d+%dj", round(realPart), round(imagPart));
    else
        printf("%d%dj", round(realPart), round(imagPart));
    end
    
    if k < N then
        printf(" , ");
    else
        printf("\n");
    end
end

K1 = 0:1:length(Xk)-1;

magnitude = abs(Xk);
subplot(3,1,2);
plot2d3(K1, magnitude);
xlabel('frequency (Hz)');
ylabel('magnitude (gain)');
title('Magnitude Spectrum');

angle = atan(imag(Xk), real(Xk));
subplot(3,1,3);
plot2d3(K1, angle);
xlabel('frequency (Hz)');
ylabel('Phase');
title('Phase Spectrum');




# OUTPUT: 
<img width="886" height="811" alt="image" src="https://github.com/user-attachments/assets/26517819-17aa-4d74-adbb-21fda361fb26" />
<img width="916" height="798" alt="image" src="https://github.com/user-attachments/assets/b1396ba6-16e0-45b0-8602-46bbf3675e2b" />



# RESULT: 
Thus, the Discrete Fourier Transform and Fast Fourier Transform of the given sequence were obtained and its magnitude and phase spectrum were plotted.

