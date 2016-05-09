function nextStates = invertPendulumDynamics(prevStateF,prevStateX1,prevStateX2)

%in radians
X1 = 1;
X2 = 2;
nextStates = zeros(1,2);

T = 0.02;
a1 = 0.0455;
a2 = 0.9091;
b1 = 0.6667;
b2 = 0.0455;
g = 9.81;

nextStates(X1) = prevStateX1 + T*prevStateX2;
nextStates(X2) = prevStateX2 + T*(g*sin(prevStateX1) + cos(prevStateX1)*(-a1*prevStateX2*prevStateX2*sin(prevStateX1) + a2*prevStateF))/(b1-b2*cos(prevStateX1)*cos(prevStateX1));

end