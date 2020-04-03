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

%%%%Matrix building
%%%%Will build 7 matrices for the 7 functions defined. 
C = [0 0 0 0 0 0 0;
   -Cap Cap 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;];

G = [1 0 0 0 0 0 0;
   -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -alpha 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+GO];
%%%%create out AC and DC voltage matrices
%%%%define our F matrix.
%V = [V1; V2; IL; V3; I3; V4; V0];
VDC=zeros(7,1);
VAC=zeros(7,1);
F=zeros(7,1);
%%%sweep from 0.1V to 10V. 10 points
%%%%%Part 3B.i
for v = 0.1:0.1:10
    F(1,1)=v;
    VDC=G\F;
    figure(1)
    plot (v, VDC(7,1),'b.') %VO
    plot(v, VDC(4,1), 'g.') %V3
    hold on
    title('DC Sweep')
    xlabel('Vin (V)')
    ylabel('Voltage (V)')   
end

%%%%%-10V to 10V sweep
%%%%Part 3B.ii
for v = -10:0.1:10
    F(1,1)=v;
    VDC=G\F;
    figure(2)
    plot (v, VDC(7,1),'b.')
    plot (v, VDC(4,1), 'g.')
    hold on
    title('DC Sweep')
    xlabel('Vin (V)')
    ylabel('Voltage (V)')   
end

%%%%%%%%% VO from AC Sweep, Gain Plot, Histogram of Gain
w = logspace(1,2,500);                  
F(1,1) = 1;
%%%%Part 3b.iii
for i = 1:length(w)
    VAC = (G+C*1j*w(i))\F;               % recalculating the voltages with AC sweeping
    figure(3)
    semilogx(w(i), abs(VAC(7,1)), 'b.')
    hold on
    xlabel('log (w)')
    ylabel('VO (V)')
    title('AC Sweep')
    gain = 20*log(abs(VAC(7,1))/F(1));    % Calculating the gain with new voltages from AC sweep
    figure(4)
    plot(i, gain, 'g.')
    hold on
    title('Gain Vo/Vin per step(dB)')
    xlabel('Step')
    ylabel('Gain (dB)')
end

% voltage gain calculation as a function of random perturbations
% on C using a normal distribution of 0.05 and w = pi
perb =  Cap + 0.05.*randn(1,1000);
w = pi;
Gain = zeros(1000,1);

for n = 1:length(Gain)
    C(1,1) = perb(n);
    C(1,2) = -perb(n);
    C(2,1) = -perb(n);
    C(2,2) = perb(n);
    VAC = (G+C*1j*w)\F;                 % Voltage calculation for AC sweeping
    Gain(n,1) = abs(VAC(7,1))/F(1);     % Gain calculation with AC voltages
end

% histogram for Gain
figure(5)
histogram(Gain,100);
title('Histogram of Gain with random perturbations')
xlabel('Gain, dB')
%%%%setting up for part 4
%%%%4a. Looking at this circuit you can see that it includes resistors,
%%%%inductors and capacitors, therefor this is an RLC circuit. 

%%%% 4b. this circuit includes inductors, capactiros and resistors. at some
%%%% frequency the inductance and capacitor reactance will equal eachother.
%%%% when this occurs there is a sharp frequency response while it was low
%%%% before. 
%Variables
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
Vin = 1;
cap = 0.25;
L = 0.2;
alpha = 100;
steps=1000;             %1000 steps
dt=1/steps;             %steps size of .03
X=1:steps;
F = [Vin;
     0;
     0;
     0;
     0;
     0;
     0];

F0 = zeros(7,1);
V0=zeros(7,steps); %%%will hold V matrix for time steps.
V1=zeros(7,1);      %%%matrix to hold V matrix while in for loop
%%%4d.ii.A
for i = 1:steps

    if i < 30
        V0(:,i) = (C./dt+G)\(F0+C*V1/dt);    
    elseif i == 30
        V0(:,i) = (C./dt+G)\(F+C*V1/dt);   
    else
        V0(:,i) = (C./dt+G)\(F+C*Vold/dt);  
    end
    Vold = V0(:, i);
    
end

figure(6)
plot(X, V0(7,:), 'b')
hold on
plot(X, V0(1,:), 'r')
title('Vin (red), VO (blue) Over Time')
xlabel('Time (ms)')
ylabel('Voltage (V)')

V2= zeros(7, steps);
freq = zeros(7,1);
V2(:,i) = (C./dt+G)\(freq+C*V1/dt);
% sin(2*pi*f*t) at a frequency of 1/30 1/ms.
Vold=zeros(7,1);
%%%%4d.ii.B
for i = 1 : steps
    freq(1) = sin(2*pi*(1/0.03)*i/steps);
    V2(:,i) = (C./dt+G)\(freq+C*Vold/dt);
    Vold = V2(:, i);  
end

figure(7)
plot(X, V2(7,:), 'b')
hold on
plot(X, V2(1,:), 'r')
title('Vin (red) Vo (blue), function sin(2pift)')
xlabel('Time (ms)')
ylabel('Voltage (V)')
%%%%4d.ii.C
%magnitude = 1, std dev. = 0.03, delay of 0.06s
V3 = zeros(7, steps);
Pulse= zeros(7,1);
V3(:,i) = (C./dt+G)\(Pulse+C*V1/dt);
Vold=zeros(7,1); 
for i = 1:steps

    Pulse(1,1) =  exp(-1/2*((i/steps-0.06)/(0.03))^2);
    V3(:,i) = (C./dt+G)\(Pulse+C*Vold/dt);   
    Vold = V3(:, i);
        
end

figure(8)
plot(X, V3(7,:), 'b')
hold on
plot(X, V3(1,:), 'r')
title('Vin (red) Vo (blue) with Guassian pulse (magnitude =1, std dev.= 30ms, delay = 60ms)')
xlabel('Time (ms)')
ylabel('Voltage (V)')

% the frequency content of the input and output signals plotted using fft() and fftshift().
%analyzing from -1/2 steps to +1/2 steps -1

% freq = (-steps/2:steps/2-1);               
% %%%%%4d.iv
% %for first signal
% fft_V1 = fft(V0.');
% ffts_V1 = fftshift(fft_V1);
% figure(9)
% plot(freq, abs(ffts_V1(:,1)), 'r')
% hold on
% plot(freq, abs(ffts_V1(:,7)), 'b')
% title('F-Domain for a) V1 (red), VO (blue)')
% xlabel('Frequency (Hz)')
% ylabel('Voltage (V)')
% 
% %Plot of Vin, Vo with second input signal in f-domain
% 
% fft_V2 = fft(V2.');
% ffts_V2 = fftshift(fft_V2);
% figure(10)
% plot(freq, abs(ffts_V2(:, 1)), 'r')
% hold on
% plot(freq, abs(ffts_V2(:, 7)), 'b')
% title('F- domain for b) V1(red) VO (blue)')
% xlabel('Frequency (Hz)')
% ylabel('Voltage (V)')
% 
% %Plot of Vin, Vo with third input signal in f-domain
% 
% fft_V3 = fft(V3.');
% ffts_V3 = fftshift(fft_V3);
% figure(11)
% plot(freq, abs(ffts_V3(:, 1)), 'r')
% hold on
% plot(freq, abs(ffts_V3(:, 7)), 'b')
% title('F- domain for c) V1 (red) VO (blue)')
% xlabel('Frequency (Hz)')
% ylabel('Voltage (V)')
% 
% freq = (-steps/2:steps/2-1);               
% %%%%%4d.V
% %for first signal
% fft_V1 = fft(V0.');
% ffts_V1 = fftshift(fft_V1);
% figure(9)
% plot(freq, abs(ffts_V1(:,1)), 'r')
% hold on
% plot(freq, abs(ffts_V1(:,7)), 'b')
% title('F-Domain for a) V1 (red), VO (blue)')
% xlabel('Frequency (Hz)')
% ylabel('Voltage (V)')
% 
% %Plot of Vin, Vo with second input signal in f-domain
% 
% fft_V2 = fft(V2.');
% ffts_V2 = fftshift(fft_V2);
% figure(10)
% plot(freq, abs(ffts_V2(:, 1)), 'r')
% hold on
% plot(freq, abs(ffts_V2(:, 7)), 'b')
% title('F- domain for b) V1(red) VO (blue)')
% xlabel('Frequency (Hz)')
% ylabel('Voltage (V)')
% 
% %Plot of Vin, Vo with third input signal in f-domain
% 
% fft_V3 = fft(V3.');
% ffts_V3 = fftshift(fft_V3);
% figure(11)
% plot(freq, abs(ffts_V3(:, 1)), 'r')
% hold on
% plot(freq, abs(ffts_V3(:, 7)), 'b')
% title('F- domain for c) V1 (red) VO (blue)')
% xlabel('Frequency (Hz)')
% ylabel('Voltage (V)')
