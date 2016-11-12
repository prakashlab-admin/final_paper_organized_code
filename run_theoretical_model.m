Nt = 2000;                 %Step Size of time
ti = 0;                   %Initial time (sec)
tf = 19.995;                  %Final time (sec)
t = linspace(ti,tf,Nt);   %Time vector (sec)

%% geometric parameters of the paperfuge system
Rr = .00025;              %Radius of rope/string (m)
Rw = .0003;               %Radius of wheel where string is attached (m)
Rfw = .040;              %Radius of the full wheel (m)
Rh = .011;               %Radius of handle (m)
l = .685/4;              %length of string (m)
w = .0003;               %thickness of wheel (m)
rho = 1000;             %density of wheel (kg/m^2)
I = (1/2)*(rho*(pi*(Rfw^2)*w))*Rfw^2;   %moment of inertia of wheel (kg m^2)
Rhw = Rw + Rh;          %constant passed into ODE

%% theoretical parameters for input force
maxForce = 40;          %maximum and minimum force applied
minForce = 0.22*maxForce;%      in the force-theta phase space (N)
minThetaStart = 35;     %thetas at which force is applied
maxThetaStart = 49;     %       in force-theta phase space (degrees)

inputForceList = [minThetaStart,minForce,maxThetaStart,maxForce];

%% free fitting parameters of the system (air, string)
aR = .003;              %air resistance coefficient from the wheel
gamma = 6;              %string twisting exponent
k = 1;                  %string effective spring constant (fixed at 1)

%% run ODE with initial conditions and plot

% x1 and v1 are initial conditions that may need to be tuned so that 
%       the ODE solution does not attenuate too quickly.
%       for example, when Rfw<.01, x1=300 works well
%                    when Rfw<.006, v1 = 7000 works well
x1 = (sind(38)*l - Rhw)/Rr; %choose x1 that corresponds to theta = 38 deg
v1 = 300;

%system output "theta" is a 2000x2 array.  
%       the first column is phi (angular position) (rad)
%       the second column is phidot (angular velocity) (rad/s)
[time,phi] = solveODE(t,Rr,Rhw,l,I,Rfw,w,aR,k,inputForceList,gamma,x1,v1);


plot(time,phi(:,2))
ylim([-2500,2500])
xlim([0,5])
ylabel('Angular Velocity (rad/s)')
xlabel('Time (s)')
title('Paperfuge Theoretical Velocity')
