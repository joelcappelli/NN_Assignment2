%%TUTE2 QU1
clear all;
close all;
clc;

%re-init
x1 = [1 1 1 1 -1 -1 1 1 1];
x2 = [-1 1 -1 -1 1 -1 -1 1 -1];
x3 = [1 1 1 -1 1 -1 -1 1 -1];

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

maxNumEpochs = 20;
L = [10 4 3];
hiddenBias = -1;

nLayers = length(L); % we'll use the number of layers often  
initWeights = cell(nLayers-1,1); % a weight matrix between each layer
 
initWeights{2} = [  0.9081 0.9485 -0.1513 -0.5108;...
                    -0.2862 -0.4820 0.6261 0.5212;...
                    0.6762 -0.4969 0.7958 0.3897...
                 ];

initWeights{1} = [  -0.0963 -0.7669 0.0939 0.1987 0.3748 -0.7527 -0.8393 0.8681 0.5382 0.3438;...
                    -0.4728 -0.5809 -0.9878 -0.1362 -0.2494 0.9467 -0.0115 -0.4997 -0.0001 0.3633;...
                    0.3310 -0.7726 -0.6481 -0.8350 0.6724 -0.9407 0.5389 -0.2807 0.4985 0.5135...
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

