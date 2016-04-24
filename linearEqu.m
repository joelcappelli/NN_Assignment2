function equ = linearEqu(point1,point2)
X=1;Y=2;
x2 = point2(X);
y2 = point2(Y);
m = linearGradient(point1,point2);
b = y2 - m*x2;
equ = [m b];
end
