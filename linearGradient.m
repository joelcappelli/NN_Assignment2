function grad = linearGradient(point1,point2)
X=1;Y=2;
x1 = point1(X);
y1 = point1(Y);
x2 = point2(X);
y2 = point2(Y);
grad = (y2 - y1)./(x2 - x1);
end