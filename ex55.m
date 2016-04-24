% ex55.m
% matlab code
% 7/4/2106
clear all;
close all
clc;

disp('(a) Max-min compositional operations ');

ra=[0.3 0.7 1.0 0.2]
ia=[0.2 0.4 0.6 0.8 1.0 0.1]
n=[0.33 0.67 1.0 0.15]
[mra,nra]=size(ra)
[mia,nia]=size(ia)
[mn,nn]=size(n)
for i=1:nra
 for j=1:nia
 p(i,j)=min(ra(i),ia(j));
 end
end
for i=1:nia
 for j=1:nn
 q(i,j)=min(ia(i),n(j));
 end
end
disp('p=raxia')
p
disp('q=iaxn')
q

% T = P o Q
% columns of Q with rows of P
% using max-min composition
T = zeros(nra,nn);

for i = 1:nra
    for j = 1:nn
        T(i,j) = max(min(p(i,:),q(:,j)'));
    end
end

T
    

ia_recovered = zeros(size(ia));
for j=1:nia
    for i = 1:nra
        ia_recovered(j) = max(ia_recovered(j),min(ra(i),p(i,j)));
    end
end 
ia_recovered
ia

ra_recovered = zeros(size(ra));
for j=1:nra
    for i = 1:nia
        ra_recovered(j) = max(ra_recovered(j),min(ia(i),p(j,i)));
    end
end
ra_recovered
ra

disp('(c) Sum-product compositional operations ');

for i=1:nra
 for j=1:nia
 p(i,j)=ra(i)*ia(j);
 end
end
disp('p=raxia')
p

ia_recovered = zeros(size(ia));
for j=1:nia
    for i = 1:nra
        ia_recovered(j) = ia_recovered(j) + ra(i)*p(i,j);
    end
end 
ia_recovered = ia_recovered/max(ia_recovered)
ia

ra_recovered = zeros(size(ra));
for j=1:nra
    for i = 1:nia
        ra_recovered(j) = ra_recovered(j) + ia(i)*p(j,i);
    end
end
ra_recovered = ra_recovered/max(ra_recovered)
ra