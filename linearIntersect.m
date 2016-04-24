function intersect = linearIntersect(line1,line2)
% line = [m, b] % gradient and y-intercept

intersect =zeros(2,1);

m1 = line1(1);
b1 = line1(2);

m2 = line2(1);
b2 = line2(2);

intersect(1) = (b2-b1)/(m1-m2);
intersect(2) = m1*intersect(1) + b1;
end