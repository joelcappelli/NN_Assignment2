%%ex2.4
clear all;
close all;
clc;

%re-init
x1 = [1;-2;0;-1];
x2 = [0;1.5;-0.5;-1];
x3 = [-1;1;0.5;-1];

setX = [x1 x2 x3];
x = setX;
% augInput = 0;
% x = [setX; augInput*ones(1,size(setX,2))];

trainingSize = size(x,2);

d1 = -1;
d2 = -1;
d3 = 1;

d = [d1 d2 d3];

eta = 0.1;

wts = [1;-1;0;0.5];

maxEpochs = 1;
patternErrors = zeros(maxEpochs,trainingSize);
wtsArray = zeros(size(wts,1),maxEpochs*trainingSize + 1);
epochs = 1;

wtsArray(:,1) = wts;

while(epochs <= maxEpochs)
    for i = 1:trainingSize
        v = wts'*x(:,i);
        z = TLU(v);    
        error = d(i)-z;
        r = error;
        patternErrors(epochs,i) = 0.5*error.*error;
        wts = wts + eta*x(:,i)*r';
        wtsArray(:,(epochs-1)*trainingSize+i + 1) = wts;
    end
    epochs = epochs +1;
end

cycleError = sum(patternErrors,2);

figure;
plot(cycleError);
title('Cycle Error');
grid on;
xlabel('Epoch');

figure;
%title('Pattern Errors');
for i = 1:trainingSize
    subplot(trainingSize,1,i);plot(patternErrors(:,i));
    grid on;
    title(strcat('Pattern ',num2str(i)));
end
xlabel('Epoch');

