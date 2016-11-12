function [ t1,y1 ] = solveODE( t,Rr,Rhw,l,I,Rfw,w,aR,k,inputForceList,gamma,x1,v1 )
%Initial conditions
initial1 = [x1, v1];

%input force parameters
minThetaStart = inputForceList(1);
minForce = inputForceList(2);
maxThetaStart = inputForceList(3);
maxForce = inputForceList(4);

%determine critical theta and phi
thetaCrit = acos(2/pi);
thetaMin = asin(Rhw/l);
phiCrit = (sin(thetaCrit)*l - Rhw)/Rr;
phiMax = sin(thetaCrit)*(l-Rhw)/Rr;

%define theoretical string force
A = k*(phiMax-phiCrit)^(gamma+1) /((gamma));
springTs = @(phi) A.*(1./(phiMax-abs(phi)).^(gamma) - 1/(phiMax)^(gamma));


[t1,y1] = ode45(@paperfugeodeSpring,t,initial1);    %ode23s also works well

function dx = paperfugeodeSpring(t,y)
%paperfugeodeSpring
%       uses current position in the paperfuge period to calculate
%       all three torques (input torque, air resistance torque, spring
%       torque) and update angular position based on ODE
    
    %use phi and phidot (y(1) and y(2)) to determine theta
    S = sign(y(1)/y(2));
    if abs(y(1))<phiCrit
        theta = asin((abs(y(1))*Rr+Rhw)/l); %general equation
    else
        %calculate theta as a function of thetaCrit based on geometry
        theta = asin(Rhw/(l - (abs(y(1))*Rr/sin(thetaCrit))));
    end
    theta = 180/pi * theta;
    
    %use certain sections of the period to determine input force
    if S>0 && theta<minThetaStart
        nowF = 0;
    elseif S>0 && theta>minThetaStart
        nowF = minForce;
    elseif S<0 && theta>maxThetaStart
        nowF = minForce;
    elseif S<0 && theta<maxThetaStart
        m = (maxForce/2)/(maxThetaStart-thetaMin);
        nowF = m*(theta-maxThetaStart) + maxForce;
    end
    
    %use input force to calculate input torque
    if abs(y(1))<phiCrit
        inputTorque = -sign(y(1))*2*Rr*nowF*(abs(y(1))*Rr + Rhw)/(l^2 - (abs(y(1))*Rr + Rhw)^2)^(1/2);
    else
        inputTorque = -sign(y(1))*2*Rr*nowF*tan(thetaCrit); %inputTorque fixed with theta of thetaCrit
    end

    %calculate torque due to air resistance
    dragWheel = aR*((4*pi/5)*Rfw^5 + w*(2*pi)*Rfw^4);
    airTorque = -sign(y(2))*dragWheel*(y(2)^2);

    %use springTs to calculate the string torque
    if S==1
        springTorque = -sign(y(1))*springTs(y(1));
    else
        springTorque = 0;
    end
    
    %update angular position and velocity using ODE
dx = [y(2);                             % Angular velocity
      (1/I)* (inputTorque + airTorque + springTorque);    % Angular acceleration
     ];
end
end