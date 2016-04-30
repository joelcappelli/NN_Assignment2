function nextStates = truckBackerUpperDynamics(prevStatePhi,prevStateTheta,prevStateX,prevStateY)

X = 1;
Y = 2;
PHI = 3;
nextStates = zeros(1,3);

b = 4;
nextStates(X) = prevStateX + cosd(prevStatePhi+prevStateTheta) + sind(prevStateTheta)*sind(prevStatePhi);
nextStates(Y) = prevStateY + sind(prevStatePhi+prevStateTheta) - sind(prevStateTheta)*cosd(prevStatePhi);
nextStates(PHI) = prevStatePhi - asind(2*sind(prevStateTheta)/b);

end
