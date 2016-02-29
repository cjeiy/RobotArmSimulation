

%% Declarations
[J,umax]=lab3robot(920607);
[J2,umax2]=lab3robot(910307);% Gives J and umax 
Lm=2; %Induction
Rm=21; %Resistance
b=1; %Friction coefficient 
J=5.5; %Moment of inertia
Ktau=38; %Material constant
Km=0.5; %Material constant
n=1/20; %Gearing factor
s=tf('s'); 

%% Assignment 1
G = ((1/(s*Lm+Rm)*Ktau*1/(J*s+b))/(1+(1/(s*Lm+Rm)*Ktau*1/(J*s+b))*Km))*1/s*n;
lab3robot(G,920607); %Checks if G is calculated correctly

disp('**Assignment 1 finished**');

%% Assignment 2
Kp=3.75532;
disp(Kp) %Prints out Kp
Gp_closedLoop=feedback(G*Kp,1); %Creates transfer function
stepinfo(Gp_closedLoop) %Use stepinfo to optimize Kp

disp('**Assignment 2 finished**');

%% Assingment 3
Gp_openLoop=G*Kp; %Open loop transferfunction including P-controller
[~,phaseMargin,~,crossOverFreq]=margin(Gp_openLoop);
disp(phaseMargin); %Prints out the phase-margin
disp(crossOverFreq) %Prints out the cross-over frequency
bandWidth = bandwidth(Gp_closedLoop);
disp(bandWidth)%Prints out the bandwidth

disp('**Assignment 3 finished**');

%% Assignment 4
figure(1)
bodeplot(G); %Plots only G without P-controller

figure(2)
bodeplot(G*Kp); %Plots G with P-controller

% Figure 1 and 2 can be compared and the questions answered!

% Interative for-loop to determine which K that makes G unstable
for K_unstable = 100:0.1:500
  G_unstable=feedback(G*K_unstable,1);
  [~,phase,~,~]=margin(G_unstable);
  if phase<=0 %If the phase<=0 we know that the system is unstable
      break
  end
end

G_unstable=feedback(G*K_unstable,1);

figure(3)
bodeplot(G_unstable); %Bode plot of the unstable system G

disp('**Assignment 4 finished**');

%% Assignment 6
figure(4)
bodeplot(Gp_openLoop); %Bode plot of the system with acceptable P-controller

[~,phase_1,~,wc]=margin(Gp_openLoop); %Get the cross-over frequency and phase-margin of the system with P-controller
disp(wc); %Prints out the cross-over frequency
disp(phase_1); %Prints out the phase-margin

wcd = 4*wc; %We want a four times faster system so we multplie the cross-over frequency with 4
disp(wcd); %Prints out the x4 cross-over frequency

figure(5)
bodeplot(G) %Bode plot of the system G where we can read the phase-margin(phase_2) and gain(Gain_wcd) at wcd  

%phase_1 = 64
%phase_2 = 25
%We want to compensate for phase_1-phase_2 = 64-25 = 39 degrees
%We want to add 6 degrees more so a total of 45 degrees to consider the lag compensator

beta = 0.175; %Figure 5.13 in the course book gives beta for 45 degrees

dB = -29.1; %The gain is read from the Bode plot of G (Figure 5) at wcd

gain_Wcd = 10^(dB/20); %The gain is calculated to be used for calculating K

K_lead_lag=1/(gain_Wcd)*sqrt(beta); %This K compensates for the -29.1dB gain margin

e1 = 0.05; %It is given that the stationary control error should be smaller than 0.05

gamma = (e1*K_lead_lag*(39.9/840)); %Gamma is calculated from the limit in the course book at page 110

tau_I = 15/wcd; %Formula from the course book

tau_D = 1/(wcd*sqrt(beta)); %Formula from the course book

%The lead and lag compensators can now be calculated as:
F_lead = ((tau_D*s)+1)/((beta*tau_D*s)+1);
F_lag = ((tau_I*s)+1)/((tau_I*s)+gamma);

%The total compensator is now:
F_compensator=K_lead_lag*F_lead*F_lag;

disp('**Assignment 6 finished**');

%% Assignment 8
Sp = 1/(1+(G*Kp)); %Sensitivity function with proportional controller
S_F = 1/(1+(G*F_compensator)); %Sensitivity function with controller F

figure(6)
bode(Sp,S_F); %Bode plot of the two sensitivity functions in the same plot

disp('**Assignment 8 finished**');

%% Assignment 9
%The robustness criterion: abs(T)<1/abs(G_delta) for all w

G_delta1 = (s+10)/40; %Modell error 1
G_delta2 = (s+10)/(4*(s+0.01)); %Model error 2

T = 1-S_F; %Complementary sensitivity function

figure(7)
bode(T,1/G_delta1); %Bode plot gives that the criterion is satisfied for G_delta1

figure(8)
bode(T,1/G_delta2); %Bode plot gives that the criterion is satisfied for w<1.34 but not for w>=1.34

disp('**Assignment 9 finished**');

%% Assignment 10
A = [0 n 0; 0 (-b/J) (Ktau/J); 0 (-Km/Lm) (-Rm/Lm)];
B = [0; 0; (1/Lm)];
C = [1 0 0];

%% Assignment 11
S=[B A*B (A*A)*B];
O=[C; C*A; C*(A*A)];

% If-else statements to know if the system is controllable and observable
if det(S)~=0
    disp('The system is controllable');
else
    disp('The system is not controllable');
end

if det(O)~=0
    disp('The system is observable');
else
    disp('The system is not observable');
end
    
%% Assignment 12
%We want one pole on the negative real axis and two poles on the "bisektris"
%The poles distance to origo gives a faster system for a relative damping i.e. poles on the "bisektris"
k = 2;
pole1 = -2;
pole2 = -k + 1*1i;
pole3 = -k - 1*1i;

Poles = [pole1 pole2 pole3];
L = place(A,B,Poles);

System = ss(A-B*L,B,C,0);

L0 = 1/dcgain(System); %s.185

figure(9)
step(L0*System);

%% Run programme
lab3robot(G,Kp,F_compensator,A,B,C,L,L0,920607);