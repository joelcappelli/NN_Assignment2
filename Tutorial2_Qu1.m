%%TUTE2 QU1
clear all;
close all;
clc;

%re-init
x1 = [1 1 1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 1];
x2 = [-1 1 -1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 1 -1 -1];
x3 = [1 1 1 1 1 -1 -1 1 1 -1 -1 1 1 1 1 1];

setX = [x1;x2;x3];
X = setX;
augInput = -1;
X = [setX augInput*ones(size(setX,1),1)];

trainingSize = size(X,2);

d1 = [1 -1 -1];
d2 = [-1 1 -1];
d3 = [-1 -1 1];

D = [d1;d2;d3];

n = 0.8;
m = 0;
smse = 1e-6;
plotFig = 1;

maxNumEpochs = 1;
L = [17 4 3];
hiddenBias = -1;

nLayers = length(L); % we'll use the number of layers often  
initWeights = cell(nLayers-1,1); % a weight matrix between each layer
 
initWeights{2} = [  -0.3124 0.6585 0.2860 -0.4668;...
                    0.0834 -0.9701 0.0923 0.9402;...
                    0.0456 -0.4516 0.8357 -0.5065...
                 ];

initWeights{1} = [  0.6880 0.6342 0.4139 0.1648 0.3914 -0.6576 -0.8402 0.5660 -0.9873 -0.8663 -0.5683 0.3181 -0.4394 0.3722 -0.7924 -0.5586 -0.6795;...
                    0.3880 -0.9558 0.4157 0.5034 -0.4410 -0.9915 -0.6738 0.4934 0.1696 0.0880 -0.2285 0.4186 0.6319 -0.3797 0.0419 -0.8894 0.4475;...
                    -0.0885 -0.6785 -0.1267 0.9831 0.5077 -0.0086 0.8066 -0.1271 0.0370 -0.2205 0.2537 0.0773 -0.2562 -0.9813 0.5518 -0.8587 -0.4060...
                  ];

% activaFn = @bipolarLog;
% activaDerivFn = @bipolarLogDeriv;
fnHandles = {@bipolarLog, @bipolarLogDeriv};

Network = backpropNN(L,n,m,smse,X,D,fnHandles,maxNumEpochs,plotFig,hiddenBias,initWeights);

% while(epochs <= maxEpochs)
%     for i = 1:trainingSize
%         v = wts'*x(:,i);
%         z = bipolarLog(v);    
%         error = d(i)-z;
%         r = error.*bipolarLog(z,'deriv');
%         patternErrors(epochs,i) = 0.5*error.*error;
%         wts = wts + eta*x(:,i)*r';
%         wtsArray(:,(epochs-1)*trainingSize+i + 1) = wts;
%     end
%     epochs = epochs +1;
% end

Network.weights{1}
Network.weights{2}

figure;
plot(Network.cycleErrors);
title('Cycle Error');
grid on;
xlabel('Epoch');

% figure;
% plot(patternErrorsOverSteps);
% title('Qu1.1 - Pattern Error Over Steps');
% grid on;
% xlabel('Step');

