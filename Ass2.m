
clear all;
close all;
clc;

%% Neural Networks Assignment 2
% Joel Cappelli
% 12137384

%% Qu1
x1 = [0.7826 0.5242];
x2 = [0.9003 -0.5377];
x3 = [-0.0871 -0.9630];
x4 = [-0.1795 0.7873];
x5 = [0.2309 0.5839];
x6 = [0.2137 -0.0280];
x7 = [-0.6475 -0.1886];
x8 = [-0.6026 0.2076];

setX = [x1;x2;x3;x4;x5;x6;x7;x8];

X = setX;
augInput = -1;
X = [setX augInput*ones(size(setX,1),1)];

d1 = 1;
d2 = -1;
d3 = 1;
d4 = -1;
d5 = 1;
d6 = -1;
d7 = 1;
d8 = -1;

D = [d1;d2;d3;d4;d5;d6;d7;d8];

n = 0.2;
m = 0;
smse = 1e-6;
plotFig = 0;

L = [3 4 1];
hiddenBias = -1;

nLayers = length(L); % we'll use the number of layers often  
initWeights = cell(nLayers-1,1); % a weight matrix between each layer
 
initWeights{2} = [ 0.3443 0.6762 -0.9607 0.3626 ];

initWeights{1} = [ -0.2410 0.4189 -0.6207;...
                   0.6636 -0.1422 -0.6131;...
                   0.0056 -0.3908 0.3644 
                 ];

% activaFn = @bipolarLog;
% activaDerivFn = @bipolarLogDeriv;
fnHandles = {@bipolarLog, @bipolarLogDeriv};

maxNumEpochs = 1;
Network = backpropNN(L,n,m,smse,X(1,:),D(1,:),fnHandles,maxNumEpochs,plotFig,hiddenBias,initWeights);

% disp('Initial weights:');
% disp(wtsArray(:,1));

fprintf('Qu1.1\n');

fprintf('Weights after 1 step (first training pattern):\n');
fprintf('Wbar(2) = \n');
disp(Network.weights{1});
fprintf('W(2) = \n');
disp(Network.weights{2});

maxNumEpochs = 500;
Network = backpropNN(L,n,m,smse,X,D,fnHandles,maxNumEpochs,plotFig,hiddenBias,initWeights);

fprintf('Qu1.2\n');

fprintf('Weights after 4000 steps (500 cycles):\n');
fprintf('Wbar(4001) = \n');
disp(Network.weights{1});
fprintf('W(4001) = \n');
disp(Network.weights{2});

figure;
plot(Network.cycleErrors);
title('Cycle Error');
grid on;
xlabel('Epoch');

%write up feedforward function 
xtest1 = [0.6263 -0.9803];
xtest2 = [0.0700 0.0500];

settestX = [xtest1;xtest2];

testX = settestX;
augInput = -1;
testX = [settestX augInput*ones(size(settestX,1),1)];

testD = [sign(prod(xtest1));sign(prod(xtest2))];
output = feedfwdMLP(testX,L,Network.weights,fnHandles,hiddenBias);

fprintf('xTV1 = \n');
disp(xtest1');
fprintf('xTV2 = \n');
disp(xtest2');

fprintf('outputTrueTV1 = \n');
disp(testD(1)');
fprintf('classifyTV1 = \n');
disp(output(1)');

fprintf('outputTrueTV2 = \n');
disp(testD(2)');
fprintf('classifyTV2 = \n');
disp(output(2)');

%% Qu 2
fprintf('\nQu2.1 Fuzzy logic Controller\n\n');
X = 1;
Y = 2;
PHI = 3;
xDesired = 10;
phiDesired = 90;

initialState = zeros(1,3);
initialState(X) = 15;
initialState(Y) = 15;
initialState(PHI) = 35;

minXPosErrorRange = 0;
maxXPosErrorRange = 20;
maxXPosError = 1;

minTruckAngleErrorRange = -90;
maxTruckAngleErrorRange = 270;
maxTruckAngleError = 1;

minControlAngleRange = -40;
maxControlAngleRange = 40;
maxOutput = 1;

initialXPosError = checkInputRange(minXPosErrorRange,maxXPosErrorRange,xDesired - initialState(X));
initialTruckAngleError = checkInputRange(minTruckAngleErrorRange,maxTruckAngleErrorRange,phiDesired - initialState(PHI));

setXPosError = {'NM','NS','ZE','PS','PM'};
truckAngleErrorTextXpos = [-60 70 90 110 240];%(minErrorRange - [-4 -0.8 0 0.8 4])./(minErrorRange-maxErrorRange);
setTruckAngleError = {'NM','NS','ZE','PS','PM'};
errorTextXpos = [1 7 10 13 19];%(minDeltaErrorRange - [-0.6 0 0.6])./(minDeltaErrorRange-maxDeltaErrorRange);
degXPosError = size(setXPosError,2);
degTruckAngleError = size(setTruckAngleError,2);

if(degXPosError > degTruckAngleError)
    setOutput = setXPosError;
else
    setOutput = setTruckAngleError;
end
contOutputTextXpos = [-35 -10 0 10 35];%(minOutputRange -[-8 -4 0 4 8])./(minOutputRange-maxOutputRange);
degOutput = size(setOutput,2);
FAM_XPosError_TruckAngleError = {'ZE','NS','NM','NM','NM';...
                                'PS','ZE','NS','NM','NM';...
                                'PM','PS','ZE','NS','NM';...
                                'PM','PM','PS','ZE','NS';...
                                'PM','PM','PM','PS','ZE'};
                            
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
logical_col = cellfun(cellfind('NS'),setXPosError);
logical_row = cellfun(cellfind('PS'),setTruckAngleError);
findOutput = cellfun(cellfind(FAM_XPosError_TruckAngleError(logical_row,logical_col)),setOutput);
test = setOutput(findOutput);

MAX = 1;
LINEAR = 2;
MIN = 3;

TruckAngleErrorFuzzySet1_Func = 'NM';
TruckAngleErrorFuzzySet1_output = [1 0 0];
TruckAngleErrorFuzzySet1_inputBreakPts = [minTruckAngleErrorRange 70 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet2_Func = 'NS';
TruckAngleErrorFuzzySet2_output = [0 1 0 0];
TruckAngleErrorFuzzySet2_inputBreakPts = [minTruckAngleErrorRange 70 90 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet3_Func = 'ZE';
TruckAngleErrorFuzzySet3_output = [0 0 1 0 0];
TruckAngleErrorFuzzySet3_inputBreakPts = [minTruckAngleErrorRange 70 90 110 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet4_Func = 'PS';
TruckAngleErrorFuzzySet4_output = [0 0 1 0];
TruckAngleErrorFuzzySet4_inputBreakPts = [minTruckAngleErrorRange 90 110 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet5_Func = 'PM';
TruckAngleErrorFuzzySet5_output = [0 0 1];
TruckAngleErrorFuzzySet5_inputBreakPts = [minTruckAngleErrorRange 110 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet_Funcs = {TruckAngleErrorFuzzySet1_Func;TruckAngleErrorFuzzySet2_Func;TruckAngleErrorFuzzySet3_Func;TruckAngleErrorFuzzySet4_Func;TruckAngleErrorFuzzySet5_Func};
TruckAngleErrorFuzzySet_output = {TruckAngleErrorFuzzySet1_output;TruckAngleErrorFuzzySet2_output;TruckAngleErrorFuzzySet3_output;TruckAngleErrorFuzzySet4_output;TruckAngleErrorFuzzySet5_output};
TruckAngleErrorFuzzySet_inputBreakPts = {TruckAngleErrorFuzzySet1_inputBreakPts;TruckAngleErrorFuzzySet2_inputBreakPts;TruckAngleErrorFuzzySet3_inputBreakPts;TruckAngleErrorFuzzySet4_inputBreakPts;TruckAngleErrorFuzzySet5_inputBreakPts};

numTruckAngleErrorFuzzySets = size(TruckAngleErrorFuzzySet_Funcs,1);

xPosErrorFuzzySet1_Func = 'NM';
xPosErrorFuzzySet1_output = [1 0 0];
xPosErrorFuzzySet1_inputBreakPts = [minXPosErrorRange 7 maxXPosErrorRange];

xPosErrorFuzzySet2_Func = 'NS';
xPosErrorFuzzySet2_output = [0 1 0 0];
xPosErrorFuzzySet2_inputBreakPts = [minXPosErrorRange 7 10 maxXPosErrorRange];

xPosErrorFuzzySet3_Func = 'ZE';
xPosErrorFuzzySet3_output = [0 0 1 0 0];
xPosErrorFuzzySet3_inputBreakPts = [minXPosErrorRange 7 10 13 maxXPosErrorRange];

xPosErrorFuzzySet4_Func = 'PS';
xPosErrorFuzzySet4_output = [0 0 1 0];
xPosErrorFuzzySet4_inputBreakPts = [minXPosErrorRange 10 13 maxXPosErrorRange];

xPosErrorFuzzySet5_Func = 'PM';
xPosErrorFuzzySet5_output = [0 0 1];
xPosErrorFuzzySet5_inputBreakPts = [minXPosErrorRange 13 maxXPosErrorRange];

xPosErrorFuzzySet_Funcs = {xPosErrorFuzzySet1_Func;xPosErrorFuzzySet2_Func;xPosErrorFuzzySet3_Func;xPosErrorFuzzySet4_Func;xPosErrorFuzzySet5_Func};
xPosErrorFuzzySet_output = {xPosErrorFuzzySet1_output;xPosErrorFuzzySet2_output;xPosErrorFuzzySet3_output;xPosErrorFuzzySet4_output;xPosErrorFuzzySet5_output};
xPosErrorFuzzySet_inputBreakPts = {xPosErrorFuzzySet1_inputBreakPts;xPosErrorFuzzySet2_inputBreakPts;xPosErrorFuzzySet3_inputBreakPts;xPosErrorFuzzySet4_inputBreakPts;xPosErrorFuzzySet5_inputBreakPts};

numxPosErrorFuzzySets = size(xPosErrorFuzzySet_Funcs,1);

contOutputFuzzySet1_Func = 'NM';
contOutputFuzzySet1_output = [1 0 0];
contOutputFuzzySet1_inputBreakPts = [minControlAngleRange -10 maxControlAngleRange];

contOutputFuzzySet2_Func = 'NS';
contOutputFuzzySet2_output = [0 1 0 0];
contOutputFuzzySet2_inputBreakPts = [minControlAngleRange -10 0 maxControlAngleRange];

contOutputFuzzySet3_Func = 'ZE';
contOutputFuzzySet3_output = [0 0 1 0 0];
contOutputFuzzySet3_inputBreakPts = [minControlAngleRange -10 0 10 maxControlAngleRange];

contOutputFuzzySet4_Func = 'PS';
contOutputFuzzySet4_output = [0 0 1 0];
contOutputFuzzySet4_inputBreakPts = [minControlAngleRange 0 10 maxControlAngleRange];

contOutputFuzzySet5_Func = 'PM';
contOutputFuzzySet5_output = [0 0 1];
contOutputFuzzySet5_inputBreakPts = [minControlAngleRange 10 maxControlAngleRange];

contOutputFuzzySet_Funcs = {contOutputFuzzySet1_Func;contOutputFuzzySet2_Func;contOutputFuzzySet3_Func;contOutputFuzzySet4_Func;contOutputFuzzySet5_Func};
contOutputFuzzySet_output = {contOutputFuzzySet1_output;contOutputFuzzySet2_output;contOutputFuzzySet3_output;contOutputFuzzySet4_output;contOutputFuzzySet5_output};
contOutputFuzzySet_inputBreakPts = {contOutputFuzzySet1_inputBreakPts;contOutputFuzzySet2_inputBreakPts;contOutputFuzzySet3_inputBreakPts;contOutputFuzzySet4_inputBreakPts;contOutputFuzzySet5_inputBreakPts};

numContOutputFuzzySets = size(contOutputFuzzySet_Funcs,1);

plotFig = 1;
if(plotFig)
    plotMembershipFunctions(contOutputFuzzySet_output,contOutputFuzzySet_inputBreakPts,setOutput,contOutputTextXpos,...
                                    xPosErrorFuzzySet_output, xPosErrorFuzzySet_inputBreakPts,setXPosError,errorTextXpos,...
                                    TruckAngleErrorFuzzySet_output, TruckAngleErrorFuzzySet_inputBreakPts,setTruckAngleError,truckAngleErrorTextXpos,...
                                    initialXPosError,initialTruckAngleError);
end
                            
truckAngleErrorFuzzyified = fuzzifiedMemFunc(TruckAngleErrorFuzzySet_Funcs,TruckAngleErrorFuzzySet_output,TruckAngleErrorFuzzySet_inputBreakPts,initialTruckAngleError,maxTruckAngleError);
xPosErrorFuzzyified = fuzzifiedMemFunc(xPosErrorFuzzySet_Funcs,xPosErrorFuzzySet_output,xPosErrorFuzzySet_inputBreakPts,initialXPosError,maxXPosError);

UMOMoutputFiringRules = determineRules(xPosErrorFuzzyified,truckAngleErrorFuzzyified,'MOM-prod',setXPosError,setTruckAngleError,setOutput,FAM_XPosError_TruckAngleError);
UCOAoutputFiringRules = determineRules(xPosErrorFuzzyified,truckAngleErrorFuzzyified,'COA-min',setXPosError,setTruckAngleError,setOutput,FAM_XPosError_TruckAngleError);
UMOM_crispVals = zeros(1,size(UMOMoutputFiringRules,2));
for i = 1:size(UMOM_crispVals,2)
    logical_row = cellfun(cellfind(UMOMoutputFiringRules{2,i}),contOutputFuzzySet_Funcs);
    output = contOutputFuzzySet_output(logical_row);
    inputBreakPts = contOutputFuzzySet_inputBreakPts(logical_row);
    inputBreakPtsArray = inputBreakPts{:};
    UMOM_crispVals(i) = mean(inputBreakPtsArray(find(output{:} ==1)));
end

disp('Fuzzification UMOM -prod');
for i = 1:size(UMOMoutputFiringRules,2)
    fprintf('Rule %d: %g/%s \n',i, UMOMoutputFiringRules{1,i},UMOMoutputFiringRules{2,i});
end

fprintf('\nMOM Method (prod) - Defuzzified Output for initialXPosError = %g, initialTruckAngleError = %g\n',initialXPosError,initialTruckAngleError);
outputAngle = dot(cell2mat(UMOMoutputFiringRules(1,:)),UMOM_crispVals)/sum(cell2mat(UMOMoutputFiringRules(1,:)));
fprintf('Output control angle (theta) = %g deg\n',outputAngle);

numSamplingPoints = 100;
% initialState = zeros(1,3);
% initialState(X) = 15;
% initialState(Y) = 15;
% initialState(PHI) = 35;
initialOutputAngle = 0;
controlAngleTheta = zeros(numSamplingPoints,1);
errorVec = zeros(numSamplingPoints,2);
stateVec = zeros(numSamplingPoints,3);
stateVec(1,:) = initialState;
controlAngleTheta(1) = initialOutputAngle;
errorVec(1,:) = [initialXPosError ...
                 initialTruckAngleError];

for i = 2:numSamplingPoints
    stateVec(i,:) = truckBackerUpperDynamics(stateVec(i-1,PHI),controlAngleTheta(i-1),stateVec(i-1,X),stateVec(i-1,Y));
    
    XPosError = checkInputRange(minXPosErrorRange,maxXPosErrorRange,xDesired - stateVec(i,X));
    TruckAngleError = checkInputRange(minTruckAngleErrorRange,maxTruckAngleErrorRange,phiDesired - stateVec(i,PHI));
    errorVec(i,:) = [XPosError TruckAngleError];
    
    truckAngleErrorFuzzyified = fuzzifiedMemFunc(TruckAngleErrorFuzzySet_Funcs,TruckAngleErrorFuzzySet_output,TruckAngleErrorFuzzySet_inputBreakPts,TruckAngleError,maxTruckAngleError);
    xPosErrorFuzzyified = fuzzifiedMemFunc(xPosErrorFuzzySet_Funcs,xPosErrorFuzzySet_output,xPosErrorFuzzySet_inputBreakPts,XPosError,maxXPosError);

    UMOMoutputFiringRules = determineRules(xPosErrorFuzzyified,truckAngleErrorFuzzyified,'MOM-prod',setXPosError,setTruckAngleError,setOutput,FAM_XPosError_TruckAngleError);
    UMOM_crispVals = zeros(1,size(UMOMoutputFiringRules,2));
    for j = 1:size(UMOM_crispVals,2)
        logical_row = cellfun(cellfind(UMOMoutputFiringRules{2,j}),contOutputFuzzySet_Funcs);
        output = contOutputFuzzySet_output(logical_row);
        inputBreakPts = contOutputFuzzySet_inputBreakPts(logical_row);
        inputBreakPtsArray = inputBreakPts{:};
        UMOM_crispVals(j) = mean(inputBreakPtsArray(find(output{:} ==1)));
    end

    controlAngleTheta(i) = dot(cell2mat(UMOMoutputFiringRules(1,:)),UMOM_crispVals)/sum(cell2mat(UMOMoutputFiringRules(1,:)));
end

figure;
subplot(2,1,1);plot(stateVec(:,X),stateVec(:,Y),'r*');
grid on;
xlabel('X');
ylabel('Y');
subplot(2,1,2);plot(controlAngleTheta,'b*');
xlabel('Sample Point');
ylabel('theta (deg)');
grid on;

figure;
subplot(2,1,1);plot(errorVec(:,X),'r*');
grid on;
xlabel('X error');
subplot(2,1,2);plot(errorVec(:,2),'b*');
xlabel('Truck Angle error');
grid on;

fprintf('\n');
disp('Fuzzification UCOA - min');
for i = 1:size(UCOAoutputFiringRules,2)
    fprintf('Rule %d: %g/%s \n',i, UCOAoutputFiringRules{1,i},UCOAoutputFiringRules{2,i});
end

UCOA_vals = cell2mat(UCOAoutputFiringRules(1,:));

UCOA_Rule1_inputBreakPts = [-10 linearPosx(linearEqu([-10 0],[0 1]),UCOA_vals(1)) linearPosx(linearEqu([0 1],[10 0]),UCOA_vals(1)) 10];
UCOA_Rule1_output = [0 UCOA_vals(1) UCOA_vals(1) 0];

figure;
for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,1,1);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,1,1);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,1,1);plot(UCOA_Rule1_inputBreakPts,UCOA_Rule1_output,'b','linewidth',3);
grid on;
title('MF1');

UCOA_Rule2_inputBreakPts = [minControlAngleRange linearPosx(linearEqu([minControlAngleRange 0],[-10 1]),UCOA_vals(2)) linearPosx(linearEqu([-10 1],[0 0]),UCOA_vals(2)) 0];
UCOA_Rule2_output = [0 UCOA_vals(2) UCOA_vals(2) 0];

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,1,2);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,1,2);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,1,2);plot(UCOA_Rule2_inputBreakPts,UCOA_Rule2_output,'b','linewidth',3);
grid on;
title('MF2');

figure;
for i = 1:numContOutputFuzzySets
str = setOutput{i};
plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
text(contOutputTextXpos(i),1,str);
hold on;
end
plot(UCOA_Rule1_inputBreakPts,UCOA_Rule1_output,'b','linewidth',3);

for i = 1:numContOutputFuzzySets
plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
hold on;
end
plot(UCOA_Rule2_inputBreakPts,UCOA_Rule2_output,'b','linewidth',3);
grid on;

% for i = 1:numContOutputFuzzySets
% plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
% hold on;
% end
% plot(UCOA_Rule3_inputBreakPts,UCOA_Rule3_output,'b','linewidth',3);
% 
% for i = 1:numContOutputFuzzySets
% plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
% hold on;
% end
% plot(UCOA_Rule4_inputBreakPts,UCOA_Rule4_output,'b','linewidth',3);
% grid on;
% 

defuzz_output = [0 UCOA_vals(2) UCOA_vals(2) UCOA_vals(1) UCOA_vals(1) 0];
defuzz_breakPts = [minControlAngleRange linearPosx(linearEqu([minControlAngleRange 0],[-10 1]),UCOA_vals(2)) linearPosx(linearEqu([-10 1],[0 0]),UCOA_vals(2)) ...
                    linearPosx(linearEqu([-10 1],[0 0]),UCOA_vals(1)) linearPosx(linearEqu([0 1],[10 0]),UCOA_vals(1)) 10];
                
plot(defuzz_breakPts,defuzz_output,'--r','linewidth',3);
grid on;
title('Defuzzification');

func1 = linearEqu([minControlAngleRange 0],[-10 1]);
func2 = [0 UCOA_vals(2)];
func3 = linearEqu([-10 1],[0 0]);
func4 = [0 UCOA_vals(1)];
func5 = linearEqu([0 1],[10 0]);

defuzz_Funcs = [func1 func2 func3 func4 func5];
areas = zeros(1,size(defuzz_breakPts,2)-1);
moments = zeros(1,size(defuzz_breakPts,2)-1);

fprintf('\n');
disp('Defuzzification Piecewise Functions');
for i = 1:(size(defuzz_breakPts,2)-1)
    if(defuzz_Funcs(2*i-1) ~= 0)
        fprintf('Function %d between  %g and %g is u(U) = %g*u + %g \n',i, defuzz_breakPts(i),defuzz_breakPts(i+1),defuzz_Funcs(2*i-1),defuzz_Funcs(2*i));
    else
        fprintf('Function %d between  %g and %g is u(U) = %g \n',i, defuzz_breakPts(i),defuzz_breakPts(i+1),defuzz_Funcs(2*i));
    end        
    areas(i) = integral(@(x)(defuzz_Funcs(2*i-1)*x + defuzz_Funcs(2*i)),defuzz_breakPts(i),defuzz_breakPts(i+1));
    moments(i) = integral(@(x)(defuzz_Funcs(2*i-1)*(x.*x) + defuzz_Funcs(2*i)*x),defuzz_breakPts(i),defuzz_breakPts(i+1));
    fprintf('Area %d = %g\n',i,areas(i));
    fprintf('Moment %d = %g\n',i,moments(i));
end

fprintf('\nCOA Method (max-min) - Defuzzified Output for initialXPosError = %g, initialTruckAngleError = %g\n',initialXPosError,initialTruckAngleError);
fprintf('Total Area = %g\n',sum(areas));
fprintf('Total Moment = %g\n',sum(moments));
outputAngle = sum(moments)/sum(areas);
fprintf('Output control angle (theta) = %g deg\n',outputAngle);

text(outputAngle,0.6,strcat('COA = ',num2str(outputAngle)));
plot([outputAngle outputAngle],[0 1],'--g','linewidth',2);
  



