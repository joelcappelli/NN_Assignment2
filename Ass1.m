
clear all;
close all;
clc;

%% Neural Networks Assignment 2
% Joel Cappelli
% 12137384

%% Qu1.1 
x1 = [0.8;0.5;0.0;0.1];
x2 = [0.2;0.1;1.3;0.9];
x3 = [0.9;0.7;0.3;0.3];
x4 = [0.2;0.7;0.8;0.2];
x5 = [1.0;0.8;0.5;0.7];
x6 = [0.0;0.2;0.3;0.6];

setX = [x1 x2 x3 x4 x5 x6];

augInput = -1;
x = [setX; augInput*ones(1,size(setX,2))];

trainingSize = size(x,2);

d1 = 1;
d2 = -1;
d3 = 1;
d4 = -1;
d5 = 1;
d6 = -1;

d = [d1 d2 d3 d4 d5 d6];

eta = 0.05;

wts = [0.0976;0.8632;-0.3296;0.3111;-0.2162];

sse = Inf;
tol = 1e-2;

maxCycles = 10;
patternErrors = zeros(maxCycles,trainingSize);
patternErrorsOverSteps = zeros(maxCycles*trainingSize,1);
wtsArray = zeros(size(wts,1),maxCycles*trainingSize + 1);
cycle = 1;

wtsArray(:,1) = wts;

while(cycle <= maxCycles)
    for i = 1:trainingSize
        v = wts'*x(:,i);
        z = TLU(v);    
        error = d(i)-z;
        r = error;
        patternErrors(cycle,i) = 0.5*error*error;
        patternErrorsOverSteps((cycle-1)*trainingSize+i) = 0.5*error*error;
        wts = wts + eta*x(:,i)*r;
        wtsArray(:,(cycle-1)*trainingSize+i+1) = wts;
    end
    cycle = cycle +1;
end

cycleError = sum(patternErrors,2);

figure;
plot(cycleError);
title('Qu1.1 - Cycle Error');
grid on;
xlabel('Cycle');

figure;
plot(patternErrorsOverSteps);
title('Qu1.1 - Pattern Error Over Steps');
grid on;
xlabel('Step');

fprintf('\nQu1.1 Discrete Perceptron Training\n');

disp('Initial weights:');
disp(wtsArray(:,1));

fprintf('Final weights after %d steps:\n',maxCycles*trainingSize);
disp(wtsArray(:,end));

figure;
for i = 1:trainingSize
    subplot(trainingSize,1,i);plot(patternErrors(:,i));
    grid on;
    title(strcat('Pattern ',num2str(i)));
end
xlabel('Cycle');

%% Qu1.2
fprintf('Qu1.2 Continuous Perceptron Training\n');

%re-init
setX = [x1 x2 x3 x4 x5 x6];
augInput = 1;
x = [setX; augInput*ones(1,size(setX,2))];

trainingSize = size(x,2);

eta = 0.5;

wts = [0.0976;0.8632;-0.3296;0.3111;-0.2162];

maxCycles = 50;
patternErrors = zeros(maxCycles,trainingSize);
patternErrorsOverSteps = zeros(maxCycles*trainingSize,1);
wtsArray = zeros(size(wts,1),maxCycles*trainingSize + 1);
cycle = 1;

wtsArray(:,1) = wts;

while(cycle <= maxCycles)
    for i = 1:trainingSize
        v = wts'*x(:,i);
        z = bipolarLog(v);    
        error = d(i)-z;
        r = error*bipolarLog(z,'deriv');
        patternErrors(cycle,i) = 0.5*error*error;
        patternErrorsOverSteps((cycle-1)*trainingSize+i) = 0.5*error*error;
        wts = wts + eta*x(:,i)*r;
        wtsArray(:,(cycle-1)*trainingSize+i+1) = wts;
    end
    cycle = cycle +1;
end

cycleError = sum(patternErrors,2);

fprintf('Desired output:\n')
disp(d);

disp('Initial weights:');
disp(wtsArray(:,1));

fprintf('Weights after %d steps (%d cycle):\n',6,1);
disp(wtsArray(:,7));

fprintf('Classification with w%d:\n',7);
classW7 = zeros(size(d));
for i = 1:trainingSize
    classW7(i) = bipolarLog(wtsArray(:,7)'*x(:,i));
end
disp(classW7);
fprintf('Cycle Error with w%d:\n',7);
disp(cycleError(1));

fprintf('Weights after %d steps (%d cycles):\n',300,50);
disp(wtsArray(:,301));

fprintf('Classification with w%d:\n',301);
classW301 = zeros(size(d));
for i = 1:trainingSize
    classW301(i) = bipolarLog(wtsArray(:,301)'*x(:,i));
end
disp(classW301);
fprintf('Cycle Error with w%d:\n',301);
disp(cycleError(50));

figure;
plot(cycleError);
title('Qu1.2 - Cycle Error');
grid on;
xlabel('Cycle');

figure;
plot(patternErrorsOverSteps);
title('Qu1.2 - Pattern Error Over Steps');
grid on;
xlabel('Step');

figure;
for i = 1:trainingSize
    subplot(trainingSize,1,i);plot(patternErrors(:,i));
    grid on;
    title(strcat('Pattern ',num2str(i)));
end
xlabel('Cycle');

%% Qu 2.1
Credit_limits = [500, 1000, 1500, 2000];
Average_account_balance = [20, 50, 100, 500, 1000, 1200];
Profits =  [-50, 0, 50, 100, 500];

[mcl,ncl]=size(Credit_limits);
[maab,naab]=size(Average_account_balance);
[mp,np]=size(Profits);

A = [
        0.8 1.0 0.6 0.2 0.0 0.0;...
        0.2 0.3 0.5 0.8 0.1 0.0;...
        0.1 0.2 0.4 0.9 0.6 0.1;...
        0.1 0.1 0.4 0.8 0.8 0.3
     ];

E = [   
        0.8 0.9 0.7 0.1 0.0;...
        0.7 1.0 0.8 0.2 0.0;...
        0.5 0.9 0.9 0.5 0.1;...
        0.2 0.5 0.7 0.8 0.9;...
        0.1 0.4 0.6 0.9 0.8;...
        0.0 0.3 0.5 0.8 0.7...
     ];

% R = A o E = (Credit_limits X Average_account_balance) o (Average_account_balance X Profits)
% columns of Q with rows of P
% using max-min composition
R_maxmincomp = zeros(ncl,np);

for i = 1:ncl
    for j = 1:np
        R_maxmincomp(i,j) = max(min(A(i,:),E(:,j)'));
    end
end

disp('Qu 2.1.1 R Fuzzy Relation - Max-min compositional');
R_maxmincomp

R_sumprodcomp = zeros(ncl,np);

for i = 1:ncl
    for j = 1:np
        R_sumprodcomp(i,j) = dot(A(i,:),E(:,j)');
    end
end

disp('Qu 2.1.2 R Fuzzy Relation - Sum-product compositional');
max_R = max(max(R_sumprodcomp));
R_sumprodcomp = R_sumprodcomp/max_R

%% Qu 2.2

fprintf('\nQu2.2 Fuzzy logic Controller\n\n');
error = 3.25;
deltaError = -0.2;

minErrorRange = -5;
maxErrorRange = 5;
maxError = 1;

minDeltaErrorRange = -2;
maxDeltaErrorRange = 2;
maxDeltaError = 1;

minOutputRange = -10;
maxOutputRange = 10;
maxOutput = 1;

setError = {'NM','NS','ZE','PS','PM'};
errorTextXpos = [-4 -0.8 0 0.8 4];%(minErrorRange - [-4 -0.8 0 0.8 4])./(minErrorRange-maxErrorRange);
setDeltaError = {'NS','ZE','PS'};
deltaErrorTextXpos = [-0.6 0 0.6];%(minDeltaErrorRange - [-0.6 0 0.6])./(minDeltaErrorRange-maxDeltaErrorRange);
degError = size(setError,2);
degDeltaError = size(setDeltaError,2);

if(degError > degDeltaError)
    setOutput = setError;
else
    setOutput = setDeltaError;
end
contOutputTextXpos = [-9 -5 0 5 9];%(minOutputRange -[-8 -4 0 4 8])./(minOutputRange-maxOutputRange);
degOutput = size(setOutput,2);
FAM_error_DeltaError = {'PM','PM','PS','ZE','NS';...
                        'PM','PS','ZE','NS','NM';...
                        'PS','ZE','NS','NM','NM'};
                   
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
logical_col = cellfun(cellfind('NS'),setError);
logical_row = cellfun(cellfind('NS'),setDeltaError);
findOutput = cellfun(cellfind(FAM_error_DeltaError(logical_row,logical_col)),setOutput);
test = setOutput(findOutput);

MAX = 1;
LINEAR = 2;
MIN = 3;

errorFuzzySet1_Func = [MAX LINEAR MIN];
errorFuzzySet1_output = [1 1 0 0];
errorFuzzySet1_inputBreakPts = [minErrorRange -4 -0.8 maxErrorRange];

errorFuzzySet2_Func = [LINEAR LINEAR MIN];
errorFuzzySet2_output = [0 1 0 0];
errorFuzzySet2_inputBreakPts = [minErrorRange -0.8 -0.2 maxErrorRange];

errorFuzzySet3_Func = [MIN LINEAR LINEAR MIN];
errorFuzzySet3_output = [0 0 1 0 0];
errorFuzzySet3_inputBreakPts = [minErrorRange -0.8 0 0.8 maxErrorRange];

errorFuzzySet4_Func = [MIN LINEAR LINEAR];
errorFuzzySet4_output = [0 0 1 0];
errorFuzzySet4_inputBreakPts = [minErrorRange 0.2 0.8 maxErrorRange];

errorFuzzySet5_Func = [MIN LINEAR MAX];
errorFuzzySet5_output = [0 0 1 1];
errorFuzzySet5_inputBreakPts = [minErrorRange 0.8 4 maxErrorRange];

errorFuzzySet_Funcs = {errorFuzzySet1_Func;errorFuzzySet2_Func;errorFuzzySet3_Func;errorFuzzySet4_Func;errorFuzzySet5_Func};
errorFuzzySet_output = {errorFuzzySet1_output;errorFuzzySet2_output;errorFuzzySet3_output;errorFuzzySet4_output;errorFuzzySet5_output};
errorFuzzySet_inputBreakPts = {errorFuzzySet1_inputBreakPts;errorFuzzySet2_inputBreakPts;errorFuzzySet3_inputBreakPts;errorFuzzySet4_inputBreakPts;errorFuzzySet5_inputBreakPts};

numErrorFuzzySets = size(errorFuzzySet_Funcs,1);

deltaErrorFuzzySet1_Func = [MAX LINEAR MIN];
deltaErrorFuzzySet1_output = [1 1 0 0];
deltaErrorFuzzySet1_inputBreakPts = [minDeltaErrorRange -0.6 0 maxDeltaErrorRange];

deltaErrorFuzzySet2_Func = [MIN LINEAR LINEAR MIN];
deltaErrorFuzzySet2_output = [0 0 1 0 0];
deltaErrorFuzzySet2_inputBreakPts = [minDeltaErrorRange -0.6 0 0.6 maxDeltaErrorRange];

deltaErrorFuzzySet3_Func = [MIN LINEAR MAX];
deltaErrorFuzzySet3_output = [0 0 1 1];
deltaErrorFuzzySet3_inputBreakPts = [minDeltaErrorRange 0 0.6 maxDeltaErrorRange];

deltaErrorFuzzySet_Funcs = {deltaErrorFuzzySet1_Func;deltaErrorFuzzySet2_Func;deltaErrorFuzzySet3_Func};
deltaErrorFuzzySet_output = {deltaErrorFuzzySet1_output;deltaErrorFuzzySet2_output;deltaErrorFuzzySet3_output};
deltaErrorFuzzySet_inputBreakPts = {deltaErrorFuzzySet1_inputBreakPts;deltaErrorFuzzySet2_inputBreakPts;deltaErrorFuzzySet3_inputBreakPts};

numdeltaErrorFuzzySets = size(deltaErrorFuzzySet_Funcs,1);

contOutputFuzzySet1_Func = [MAX LINEAR MIN];
contOutputFuzzySet1_output = [1 1 0 0];
contOutputFuzzySet1_inputBreakPts = [minOutputRange -9 -5 maxOutputRange];

contOutputFuzzySet2_Func = [LINEAR LINEAR MIN];
contOutputFuzzySet2_output = [0 1 0 0];
contOutputFuzzySet2_inputBreakPts = [minOutputRange -5 -1 maxOutputRange];

contOutputFuzzySet3_Func = [MIN LINEAR LINEAR MIN];
contOutputFuzzySet3_output = [0 0 1 0 0];
contOutputFuzzySet3_inputBreakPts = [minOutputRange -5 0 5 maxOutputRange];

contOutputFuzzySet4_Func = [MIN LINEAR LINEAR];
contOutputFuzzySet4_output = [0 0 1 0];
contOutputFuzzySet4_inputBreakPts = [minOutputRange 1 5 maxOutputRange];

contOutputFuzzySet5_Func = [MIN LINEAR MAX];
contOutputFuzzySet5_output = [0 0 1 1];
contOutputFuzzySet5_inputBreakPts = [minOutputRange 5 9 maxOutputRange];

contOutputFuzzySet_Funcs = {contOutputFuzzySet1_Func;contOutputFuzzySet2_Func;contOutputFuzzySet3_Func;contOutputFuzzySet4_Func;contOutputFuzzySet5_Func};
contOutputFuzzySet_output = {contOutputFuzzySet1_output;contOutputFuzzySet2_output;contOutputFuzzySet3_output;contOutputFuzzySet4_output;contOutputFuzzySet5_output};
contOutputFuzzySet_inputBreakPts = {contOutputFuzzySet1_inputBreakPts;contOutputFuzzySet2_inputBreakPts;contOutputFuzzySet3_inputBreakPts;contOutputFuzzySet4_inputBreakPts;contOutputFuzzySet5_inputBreakPts};

numContOutputFuzzySets = size(contOutputFuzzySet_Funcs,1);

figure;
for i = 1:numErrorFuzzySets
str = setError{i};
subplot(3,1,1); plot(errorFuzzySet_inputBreakPts{i},errorFuzzySet_output{i},'b');
subplot(3,1,1); text(errorTextXpos(i),1,str);
hold on;
end
subplot(3,1,1);plot([error error],[0 1],'--r','linewidth',2);
grid on;
ylabel('Error');

for i = 1:numdeltaErrorFuzzySets
str = setDeltaError{i};
subplot(3,1,2); plot(deltaErrorFuzzySet_inputBreakPts{i},deltaErrorFuzzySet_output{i},'b');
subplot(3,1,2);text(deltaErrorTextXpos(i),1,str);
hold on;
end
subplot(3,1,2);plot([deltaError deltaError],[0 1],'--k','linewidth',2);
grid on;
ylabel('Delta Error');

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(3,1,3); plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(3,1,3);text(contOutputTextXpos(i),1,str);
hold on;
end
grid on;
ylabel('Controller Output');

inputs = minErrorRange:0.01:maxErrorRange;
output = zeros(size(inputs));
numinputs = size(inputs,2);

fuzzySet_Error_PM = polyval(linearEqu([0.8 0],[4 1]),error);
fuzzySet_Error_PS = polyval(linearEqu([0.8 1],[5 0]),error);

fuzzySet_deltaError_NS = polyval(linearEqu([-0.6 1],[0 0]),deltaError);
fuzzySet_deltaError_ZE = polyval(linearEqu([-0.6 0],[0 1]),deltaError);

logical_col = cellfun(cellfind('PS'),setError);
logical_row = cellfun(cellfind('NS'),setDeltaError);
findOutput = cellfun(cellfind(FAM_error_DeltaError(logical_row,logical_col)),setOutput);
Rule1_output = setOutput(findOutput);
logical_col = cellfun(cellfind('PS'),setError);
logical_row = cellfun(cellfind('ZE'),setDeltaError);
findOutput = cellfun(cellfind(FAM_error_DeltaError(logical_row,logical_col)),setOutput);
Rule2_output = setOutput(findOutput);
logical_col = cellfun(cellfind('PM'),setError);
logical_row = cellfun(cellfind('NS'),setDeltaError);
findOutput = cellfun(cellfind(FAM_error_DeltaError(logical_row,logical_col)),setOutput);
Rule3_output = setOutput(findOutput);
logical_col = cellfun(cellfind('PM'),setError);
logical_row = cellfun(cellfind('ZE'),setDeltaError);
findOutput = cellfun(cellfind(FAM_error_DeltaError(logical_row,logical_col)),setOutput);
Rule4_output = setOutput(findOutput);

Rule_output = {Rule1_output{:},Rule2_output{:},Rule3_output{:},Rule4_output{:}};

UMOM_Rule1_E_PS_deltaE_NS = fuzzySet_Error_PS*fuzzySet_deltaError_NS;
UMOM_Rule2_E_PS_deltaE_ZE = fuzzySet_Error_PS*fuzzySet_deltaError_ZE;
UMOM_Rule3_E_PM_deltaE_NS = fuzzySet_Error_PM*fuzzySet_deltaError_NS;
UMOM_Rule4_E_PM_deltaE_ZE = fuzzySet_Error_PM*fuzzySet_deltaError_ZE;

UMOM_vals = [UMOM_Rule1_E_PS_deltaE_NS,UMOM_Rule2_E_PS_deltaE_ZE,UMOM_Rule3_E_PM_deltaE_NS,UMOM_Rule4_E_PM_deltaE_ZE];

UCOA_Rule1_E_PS_deltaE_NS = min(fuzzySet_Error_PS,fuzzySet_deltaError_NS);
UCOA_Rule2_E_PS_deltaE_ZE = min(fuzzySet_Error_PS,fuzzySet_deltaError_ZE);
UCOA_Rule3_E_PM_deltaE_NS = min(fuzzySet_Error_PM,fuzzySet_deltaError_NS);
UCOA_Rule4_E_PM_deltaE_ZE = min(fuzzySet_Error_PM,fuzzySet_deltaError_ZE);

UCOA_vals = [UCOA_Rule1_E_PS_deltaE_NS,UCOA_Rule2_E_PS_deltaE_ZE,UCOA_Rule3_E_PM_deltaE_NS,UCOA_Rule4_E_PM_deltaE_ZE];

disp('Fuzzification UMOM - sum-prod');
for i = 1:size(UMOM_vals,2)
    fprintf('Rule %d: %g/%s \n',i, UMOM_vals(i),Rule_output{i});
end

fprintf('\nMOM Method (sum-prod) - Defuzzified Output for Error = %g, deltaError = %g\n',error,deltaError);
outputVoltage = dot(UMOM_vals,[0 -5 -5 -9])/sum(UMOM_vals);
fprintf('Output voltage = %gV\n',outputVoltage);


fprintf('\n');
disp('Fuzzification UCOA - max-min');
for i = 1:size(UMOM_vals,2)
    fprintf('Rule %d: %g/%s \n',i, UCOA_vals(i),Rule_output{i});
end

UCOA_Rule1_inputBreakPts = [-5 linearPosx(linearEqu([-5 0],[0 1]),UCOA_vals(1)) linearPosx(linearEqu([5 0],[0 1]),UCOA_vals(1)) 5];
UCOA_Rule1_output = [0 UCOA_vals(1) UCOA_vals(1) 0];

figure;
for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,2,1);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,2,1);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,2,1);plot(UCOA_Rule1_inputBreakPts,UCOA_Rule1_output,'b','linewidth',3);
grid on;
title('MF1');

UCOA_Rule2_inputBreakPts = [-10 linearPosx(linearEqu([-10 0],[-5 1]),UCOA_vals(2)) linearPosx(linearEqu([-5 1],[-1 0]),UCOA_vals(2)) -1];
UCOA_Rule2_output = [0 UCOA_vals(2) UCOA_vals(2) 0];

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,2,2);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,2,2);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,2,2);plot(UCOA_Rule2_inputBreakPts,UCOA_Rule2_output,'b','linewidth',3);
grid on;
title('MF2');

UCOA_Rule3_inputBreakPts = [-10 linearPosx(linearEqu([-10 0],[-5 1]),UCOA_vals(3)) linearPosx(linearEqu([-5 1],[-1 0]),UCOA_vals(3)) -1];
UCOA_Rule3_output = [0 UCOA_vals(3) UCOA_vals(3) 0];

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,2,3);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,2,3);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,2,3);plot(UCOA_Rule3_inputBreakPts,UCOA_Rule3_output,'b','linewidth',3);
grid on;
title('MF3');

UCOA_Rule4_inputBreakPts = [-10 linearPosx(linearEqu([-9 1],[-5 0]),UCOA_vals(4)) -5];
UCOA_Rule4_output = [UCOA_vals(4) UCOA_vals(4) 0];

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,2,4);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,2,4);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,2,4);plot(UCOA_Rule4_inputBreakPts,UCOA_Rule4_output,'b','linewidth',3);
grid on;
title('MF4');

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

for i = 1:numContOutputFuzzySets
plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
hold on;
end
plot(UCOA_Rule3_inputBreakPts,UCOA_Rule3_output,'b','linewidth',3);

for i = 1:numContOutputFuzzySets
plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
hold on;
end
plot(UCOA_Rule4_inputBreakPts,UCOA_Rule4_output,'b','linewidth',3);
grid on;

defuzz_output = [UCOA_vals(4) UCOA_vals(4) UCOA_vals(2) UCOA_vals(2) UCOA_vals(1) UCOA_vals(1) 0];
defuzz_breakPts = [-10 linearPosx(linearEqu([-9 1],[-5 0]),UCOA_vals(4)) linearPosx(linearEqu([-9 1],[-5 0]),UCOA_vals(2)) ...
                    linearPosx(linearEqu([-5 1],[-1 0]),UCOA_vals(2)) linearPosx(linearEqu([-5 1],[-1 0]),UCOA_vals(1)) linearPosx(linearEqu([0 1],[5 0]),UCOA_vals(1)) 5];
                
plot(defuzz_breakPts,defuzz_output,'--r','linewidth',3);
grid on;
title('Defuzzification');

func1 = [0 UCOA_vals(4)];
func2 = linearEqu([-9 1],[-5 0]);
func3 = [0 UCOA_vals(2)];
func4 = linearEqu([-5 1],[-1 0]);
func5 = [0 UCOA_vals(1)];
func6 = linearEqu([0 1],[5 0]);

defuzz_Funcs = [func1 func2 func3 func4 func5 func6];
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

fprintf('\nCOA Method (max-min) - Defuzzified Output for Error = %g, deltaError = %g\n',error,deltaError);
fprintf('Total Area = %g\n',sum(areas));
fprintf('Total Moment = %g\n',sum(moments));
outputVoltage = sum(moments)/sum(areas);
fprintf('Output voltage = %gV\n',outputVoltage);
text(outputVoltage,0.6,strcat('COA = ',num2str(outputVoltage)));
plot([outputVoltage outputVoltage],[0 1],'--g','linewidth',2);

    



