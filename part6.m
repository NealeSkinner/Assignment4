%%%V = ?I3 + ?I3^2 + ?I3^3
%%% to properly model this we would need to add values for betta and gamma.
%%% in the G matrix where alpha once was, we would include betta and gamma
%%% as well. i would then be able to go through the same format that was
%%% seen in part 1 to solve the circuit. 

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
beta=100e5;
gamma=100e8;
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
    0 0 0 0 0 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+GO];
I3=0;
J = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 (-alpha - beta*I3.^2 - gamma*I3.^3);];
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
