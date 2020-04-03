%%%%Variable Definitions
R1 = 1;
G1 = 1/R1;
R2 = 2;
G2 = 1/R2;
R3 = 10;
G3 = 1/R3;
R4 = 0.1;
G4 = 1/R4;
RO = 1000;
GO = 1/RO;
L = 0.2;
Cap = 0.25;
alpha = 100;
Cn = 0.00001;
Cn2= Cn/5.2;
Cn3= Cn*5.2;            %%%seems to have a drastic change just after *5
Vin=1;
steps=1000;             %1000 steps
steps2=steps/5;
steps3=steps*5;
dt1=1/steps;             %steps size of .03
dt2=1/steps2;
dt3=1/steps3;
X=1:steps;
%%%%Matrix building
%%%%Will build 7x7 matrices for the 7 Variables defined. 
%%%%the G and C matrices are different now for including noise figures and
%%%%currents. 
%[V1; V2; IL; V3; I3; V4; V0];
C1 = [0 0 0 0 0 0 0;
    -Cap Cap 0 0 0 0 0;
     0 0 -L 0 0 0 0;
     0 0 0 -Cn 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 -Cn 0 0 0;
     0 0 0 0 0 0 0;];

C2 = [0 0 0 0 0 0 0;
    -Cap Cap 0 0 0 0 0;
     0 0 -L 0 0 0 0;
     0 0 0 -Cn2 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 -Cn2 0 0 0;
     0 0 0 0 0 0 0;];

C3 = [0 0 0 0 0 0 0;
    -Cap Cap 0 0 0 0 0;
     0 0 -L 0 0 0 0;
     0 0 0 -Cn3 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 -Cn3 0 0 0;
     0 0 0 0 0 0 0;];

G = [1 0 0 0 0 0 0;
    -G2 G1+G2 -1 0 0 0 0;
      0 1 0 -1 0 0 0;
      0 0 -1 G3 0 0 0;
      0 0 0 0 -alpha 1 0;
      0 0 0 G3 -1 0 0;
      0 0 0 0 0 -G4 G4+GO];
F = [Vin;
     0;
     0;
     0;
     0;
     0;
     0];
V1=zeros(7,1);      %%%matrix to hold V matrix while in for loop
Pulse= zeros(7,1);
Vold=zeros(7,1);
V1(:,1) = (C1./dt1+G)\(Pulse+C1*Vold/dt1);
for i = 2:steps
    Pulse(4,1) = 0.001*randn();     %%%random noise function for In
    Pulse(1,1) =  exp(-1/2*((i/steps-0.06)/(0.03))^2);      %gaussian pulse
    V1(:,i) = (C1./dt1+G)\(Pulse+C1*Vold/dt1);   
    Vold = V1(:, i); 
end
figure(1)
plot(X, V1(7,:), 'b')
hold on
plot(X, V1(1,:), 'r')
%Fourier Transform plot
freq = (-steps/2:steps/2-1);               
fft_V1 = fft(V1.');
ffts_V1 = fftshift(fft_V1);
figure(2)
plot(freq, abs(ffts_V1(:, 1)), 'b')
hold on
plot(freq, abs(ffts_V1(:, 7)), 'r')
title('Fourier-Transform Plot of VO')
xlabel('frequency (Hz)')
ylabel('Voltage (V)')
grid on
%%%using different Cn values
V2=zeros(7,1);      %%%matrix to hold V matrix while in for loop
Pulse2= zeros(7,1);
Vold=zeros(7,1);
V2(:,1) = (C2./dt1+G)\(Pulse2+C2*Vold/dt1);
for i = 2:steps
    Pulse2(4,1) = 0.001*randn();     %%%random noise function for In
    Pulse2(1,1) =  exp(-1/2*((i/steps-0.06)/(0.03))^2);      %gaussian pulse
    V2(:,i) = (C2./dt1+G)\(Pulse2+C2*Vold/dt1);   
    Vold = V2(:, i); 
end
figure(3)
plot(X, V2(7,:), 'b')
hold on
plot(X, V2(1,:), 'r')
title('Vin (red) VO (blue), Cn Smaller than Original')
xlabel('Time (ms)')
ylabel('Voltage (V)')

V3=zeros(7,1);      %%%matrix to hold V matrix while in for loop
Pulse= zeros(7,1);
Vold=zeros(7,1);
V3(:,1) = (C3./dt1+G)\(Pulse+C3*Vold/dt1);
for i = 2:steps
    Pulse(4,1) = 0.001*randn();     %%%random noise function for In
    Pulse(1,1) =  exp(-1/2*((i/steps-0.06)/(0.03))^2);      %gaussian pulse
    V3(:,i) = (C3./dt1+G)\(Pulse+C3*Vold/dt1);   
    Vold = V3(:, i); 
end
figure(4)
plot(X, V3(7,:), 'b')
hold on
plot(X, V3(1,:), 'r')
title('Vin (red) VO (blue), Cn larger than Original')
xlabel('Time (ms)')
ylabel('Voltage (V)')

%%%%Plots for changing time steps
V4=zeros(7,1);      %%%matrix to hold V matrix while in for loop
Pulse= zeros(7,1);
Vold=zeros(7,1);
V4(:,1) = (C1./dt2+G)\(Pulse+C1*Vold/dt2);
for i = 2:steps2
    Pulse(4,1) = 0.001*randn();     %%%random noise function for In
    Pulse(1,1) =  exp(-1/2*((i/steps2-0.06)/(0.03))^2);      %gaussian pulse
    V4(:,i) = (C1./dt2+G)\(Pulse+C1*Vold/dt2);   
    Vold = V4(:, i); 
end
X2=1:steps2;
figure(5)
plot(X2, V4(7,:), 'b')
hold on
plot(X2, V4(1,:), 'r')
title('Vin (red) VO (blue), time step larger than oringial')
xlabel('Time (ms)')
ylabel('Voltage (V)')


V5=zeros(7,1);      %%%matrix to hold V matrix while in for loop
Pulse= zeros(7,1);
Vold=zeros(7,1);
V5(:,1) = (C1./dt3+G)\(Pulse+C1*Vold/dt3);
for i = 2:steps3
    Pulse(4,1) = 0.001*randn();     %%%random noise function for In
    Pulse(1,1) =  exp(-1/2*((i/steps3-0.06)/(0.03))^2);      %gaussian pulse
    V5(:,i) = (C1./dt3+G)\(Pulse+C1*Vold/dt3);   
    Vold = V5(:, i); 
end
X2=1:steps3;
figure(6)
plot(X2, V5(7,:), 'b')
hold on
plot(X2, V5(1,:), 'r')
title('Vin (red) VO (blue), time step smaller than orignal')
xlabel('Time (ms)')
ylabel('Voltage (V)')
