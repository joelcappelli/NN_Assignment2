
clear all
close all;
clc
%% Qu 2
fprintf('\nTutorial Qu2 Fuzzy logic Controller\n\n');

X1 = 1;
X2 = 2;
xDesired = 0;
phiDesired = 0;

initialState = zeros(1,2);
initialState(X1) = -0.12;
initialState(X2) = -0.3;

minXPosErrorRange = -0.2;
maxXPosErrorRange = 0.2;
maxXPosError = 1;

minTruckAngleErrorRange = -1;
maxTruckAngleErrorRange = 1;
maxTruckAngleError = 1;

minControlAngleRange = -10;
maxControlAngleRange = 10;
maxOutput = 1;

initialXPosError = initialState(X1);%checkInputRange(minXPosErrorRange,maxXPosErrorRange,xDesired - initialState(X1));
initialTruckAngleError = initialState(X2);%checkInputRange(minTruckAngleErrorRange,maxTruckAngleErrorRange,phiDesired - initialState(PHI));

setXPosError = {'NM','NS','ZE','PS','PM'};
truckAngleErrorTextXpos = [-0.9 -0.5 0 0.5 0.9];%(minErrorRange - [-4 -0.8 0 0.8 4])./(minErrorRange-maxErrorRange);
setTruckAngleError = {'NM','NS','ZE','PS','PM'};
errorTextXpos = [-0.15 -0.1 0 0.1 0.15];%(minDeltaErrorRange - [-0.6 0 0.6])./(minDeltaErrorRange-maxDeltaErrorRange);
degXPosError = size(setXPosError,2);
degTruckAngleError = size(setTruckAngleError,2);

if(degXPosError > degTruckAngleError)
    setOutput = setXPosError;
else
    setOutput = setTruckAngleError;
end
contOutputTextXpos = [-9 -5 0 5 9];%(minOutputRange -[-8 -4 0 4 8])./(minOutputRange-maxOutputRange);
degOutput = size(setOutput,2);
FAM_XPosError_TruckAngleError = {'PM','PM','PM','PS','ZE';...
                                'PM','PM','PS','ZE','NS';...
                                'PM','PS','ZE','NS','NM';...
                                'PS','ZE','NS','NM','NM';...
                                'ZE','NS','NM','NM','NM'};
                            
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
TruckAngleErrorFuzzySet1_inputBreakPts = [minTruckAngleErrorRange -0.5 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet2_Func = 'NS';
TruckAngleErrorFuzzySet2_output = [0 1 0 0];
TruckAngleErrorFuzzySet2_inputBreakPts = [minTruckAngleErrorRange -0.5 0 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet3_Func = 'ZE';
TruckAngleErrorFuzzySet3_output = [0 0 1 0 0];
TruckAngleErrorFuzzySet3_inputBreakPts = [minTruckAngleErrorRange -0.5 0 0.5 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet4_Func = 'PS';
TruckAngleErrorFuzzySet4_output = [0 0 1 0];
TruckAngleErrorFuzzySet4_inputBreakPts = [minTruckAngleErrorRange 0 0.5 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet5_Func = 'PM';
TruckAngleErrorFuzzySet5_output = [0 0 1];
TruckAngleErrorFuzzySet5_inputBreakPts = [minTruckAngleErrorRange 0.5 maxTruckAngleErrorRange];

TruckAngleErrorFuzzySet_Funcs = {TruckAngleErrorFuzzySet1_Func;TruckAngleErrorFuzzySet2_Func;TruckAngleErrorFuzzySet3_Func;TruckAngleErrorFuzzySet4_Func;TruckAngleErrorFuzzySet5_Func};
TruckAngleErrorFuzzySet_output = {TruckAngleErrorFuzzySet1_output;TruckAngleErrorFuzzySet2_output;TruckAngleErrorFuzzySet3_output;TruckAngleErrorFuzzySet4_output;TruckAngleErrorFuzzySet5_output};
TruckAngleErrorFuzzySet_inputBreakPts = {TruckAngleErrorFuzzySet1_inputBreakPts;TruckAngleErrorFuzzySet2_inputBreakPts;TruckAngleErrorFuzzySet3_inputBreakPts;TruckAngleErrorFuzzySet4_inputBreakPts;TruckAngleErrorFuzzySet5_inputBreakPts};

numTruckAngleErrorFuzzySets = size(TruckAngleErrorFuzzySet_Funcs,1);

xPosErrorFuzzySet1_Func = 'NM';
xPosErrorFuzzySet1_output = [1 0 0];
xPosErrorFuzzySet1_inputBreakPts = [minXPosErrorRange -0.1 maxXPosErrorRange];

xPosErrorFuzzySet2_Func = 'NS';
xPosErrorFuzzySet2_output = [0 1 0 0];
xPosErrorFuzzySet2_inputBreakPts = [minXPosErrorRange -0.1 0 maxXPosErrorRange];

xPosErrorFuzzySet3_Func = 'ZE';
xPosErrorFuzzySet3_output = [0 0 1 0 0];
xPosErrorFuzzySet3_inputBreakPts = [minXPosErrorRange -0.1 0 0.1 maxXPosErrorRange];

xPosErrorFuzzySet4_Func = 'PS';
xPosErrorFuzzySet4_output = [0 0 1 0];
xPosErrorFuzzySet4_inputBreakPts = [minXPosErrorRange 0 0.1 maxXPosErrorRange];

xPosErrorFuzzySet5_Func = 'PM';
xPosErrorFuzzySet5_output = [0 0 1];
xPosErrorFuzzySet5_inputBreakPts = [minXPosErrorRange 0.1 maxXPosErrorRange];

xPosErrorFuzzySet_Funcs = {xPosErrorFuzzySet1_Func;xPosErrorFuzzySet2_Func;xPosErrorFuzzySet3_Func;xPosErrorFuzzySet4_Func;xPosErrorFuzzySet5_Func};
xPosErrorFuzzySet_output = {xPosErrorFuzzySet1_output;xPosErrorFuzzySet2_output;xPosErrorFuzzySet3_output;xPosErrorFuzzySet4_output;xPosErrorFuzzySet5_output};
xPosErrorFuzzySet_inputBreakPts = {xPosErrorFuzzySet1_inputBreakPts;xPosErrorFuzzySet2_inputBreakPts;xPosErrorFuzzySet3_inputBreakPts;xPosErrorFuzzySet4_inputBreakPts;xPosErrorFuzzySet5_inputBreakPts};

numxPosErrorFuzzySets = size(xPosErrorFuzzySet_Funcs,1);

contOutputFuzzySet1_Func = 'NM';
contOutputFuzzySet1_output = [1 0 0];
contOutputFuzzySet1_inputBreakPts = [minControlAngleRange -5 maxControlAngleRange];

contOutputFuzzySet2_Func = 'NS';
contOutputFuzzySet2_output = [0 1 0 0];
contOutputFuzzySet2_inputBreakPts = [minControlAngleRange -5 0 maxControlAngleRange];

contOutputFuzzySet3_Func = 'ZE';
contOutputFuzzySet3_output = [0 0 1 0 0];
contOutputFuzzySet3_inputBreakPts = [minControlAngleRange -5 0 5 maxControlAngleRange];

contOutputFuzzySet4_Func = 'PS';
contOutputFuzzySet4_output = [0 0 1 0];
contOutputFuzzySet4_inputBreakPts = [minControlAngleRange 0 5 maxControlAngleRange];

contOutputFuzzySet5_Func = 'PM';
contOutputFuzzySet5_output = [0 0 1];
contOutputFuzzySet5_inputBreakPts = [minControlAngleRange 5 maxControlAngleRange];

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

fprintf('Qu2.1 COA Strategy\n\n');
UCOAoutputFiringRules = determineRules(xPosErrorFuzzyified,truckAngleErrorFuzzyified,'COA-min',setXPosError,setTruckAngleError,setOutput,FAM_XPosError_TruckAngleError);

disp('Fuzzification UCOA - min');
for i = 1:size(UCOAoutputFiringRules,2)
    fprintf('Rule %d: %g/%s \n',i, UCOAoutputFiringRules{1,i},UCOAoutputFiringRules{2,i});
end

UCOA_vals = cell2mat(UCOAoutputFiringRules(1,:));

UCOA_Rule1_inputBreakPts = [minControlAngleRange 5 linearPosx(linearEqu([5 0],[maxControlAngleRange 1]),UCOA_vals(1)) maxControlAngleRange];
UCOA_Rule1_output = [0 0 UCOA_vals(1) UCOA_vals(1)];

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

UCOA_Rule2_inputBreakPts = [minControlAngleRange 5 linearPosx(linearEqu([5 0],[maxControlAngleRange 1]),UCOA_vals(2)) maxControlAngleRange];
UCOA_Rule2_output = [0 0 UCOA_vals(2) UCOA_vals(2)];

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,2,2);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,2,2);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,2,2);plot(UCOA_Rule2_inputBreakPts,UCOA_Rule2_output,'b','linewidth',3);
grid on;
title('MF2');

UCOA_Rule3_inputBreakPts = [minControlAngleRange 5 linearPosx(linearEqu([5 0],[maxControlAngleRange 1]),UCOA_vals(3)) maxControlAngleRange];
UCOA_Rule3_output = [0 0 UCOA_vals(3) UCOA_vals(3)];

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(2,2,3);plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(2,2,3);text(contOutputTextXpos(i),1,str);
hold on;
end
subplot(2,2,3);plot(UCOA_Rule3_inputBreakPts,UCOA_Rule3_output,'b','linewidth',3);
grid on;
title('MF3');

UCOA_Rule4_inputBreakPts = [minControlAngleRange 0 linearPosx(linearEqu([0 0],[5 1]),UCOA_vals(4)) linearPosx(linearEqu([5 1],[maxControlAngleRange 0]),UCOA_vals(4)) maxControlAngleRange];
UCOA_Rule4_output = [0 0 UCOA_vals(4) UCOA_vals(4) 0];

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
grid on;

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

defuzz_output = [0 UCOA_vals(4) UCOA_vals(4) UCOA_vals(3) UCOA_vals(3)];
defuzz_breakPts = [0 linearPosx(linearEqu([0 0],[5 1]),UCOA_vals(4)) linearPosx(linearEqu([5 0],[maxControlAngleRange 1]),UCOA_vals(4))...
    linearPosx(linearEqu([5 0],[maxControlAngleRange 1]),UCOA_vals(3)) maxControlAngleRange];

plot(defuzz_breakPts,defuzz_output,'--r','linewidth',3);
grid on;
title('Defuzzification');

% y = mx + b 
% [m b]
func1 = linearEqu([0 0],[5 1]);
func2 = [0 UCOA_vals(4)];
func3 = linearEqu([5 0],[maxControlAngleRange 1]);
func4 = [0 UCOA_vals(3)];
% func5 = linearEqu([0 1],[10 0]);

defuzz_Funcs = [func1 func2 func3 func4];
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

fprintf('\nCOA Method (max-min) - Defuzzified Output for initialXPosError = %g m, initialTruckAngleError = %g deg\n',initialXPosError,initialTruckAngleError);
fprintf('Total Area = %g\n',sum(areas));
fprintf('Total Moment = %g\n',sum(moments));
outputAngle = sum(moments)/sum(areas);
fprintf('theta(%d) = %g deg\n',1,outputAngle);
state2Vec_COA = invertPendulumDynamics(outputAngle,initialState(X1),initialState(X2));
fprintf('[x1(%d), x2(%d)] = [%g, %g]\n\n',2,2,state2Vec_COA(1,X1),state2Vec_COA(1,X2));

text(outputAngle,0.6,strcat('COA = ',num2str(outputAngle)));
plot([outputAngle outputAngle],[0 1],'--g','linewidth',2);
  
fprintf('\nQu2.2 MOM Strategy\n\n');

numSamplingPoints = 100;
% initialState = zeros(1,3);
% initialState(X1) = 15;
% initialState(X2) = 15;
% initialState(PHI) = 35;
initialOutputAngle = 0;
controlAngleTheta = zeros(numSamplingPoints,1);
errorVec = zeros(numSamplingPoints,2);
stateVec = zeros(numSamplingPoints,2);
stateVec(1,:) = initialState;

for i = 1:(numSamplingPoints-1)
    
    XPosError = checkInputRange(minXPosErrorRange,maxXPosErrorRange,stateVec(i,X1));
    TruckAngleError = checkInputRange(minTruckAngleErrorRange,maxTruckAngleErrorRange,stateVec(i,X2));
    errorVec(i,:) = [xDesired - XPosError,phiDesired - TruckAngleError];
    
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
    stateVec(i+1,:) = invertPendulumDynamics(controlAngleTheta(i),stateVec(i,X1),stateVec(i,X2));

    if(i==1 || i ==2)
        disp(strcat('Fuzzification UMOM -prod: State: ',num2str(i)));
        for j = 1:size(UMOMoutputFiringRules,2)
            fprintf('Rule %d: %g/%s \n',j, UMOMoutputFiringRules{1,j},UMOMoutputFiringRules{2,j});
        end

        fprintf('MOM Method (prod) - Defuzzified Output for xPosError = %g m, truckAngleError = %g deg\n',XPosError,TruckAngleError);
        fprintf('theta(%d) = %g deg\n',i,controlAngleTheta(i));
        fprintf('[x1(%d), x2(%d)] = [%g, %g]\n\n',i+1,i+1,stateVec(i+1,X1),stateVec(i+1,X2));
    end

end

figure;
subplot(2,1,1);plot(stateVec(:,X1),stateVec(:,X2),'r*');
grid on;
xlabel('X1');
ylabel('X2');
subplot(2,1,2);plot(controlAngleTheta,'b*');
xlabel('Sample Point');
ylabel('theta (deg)');
grid on;

figure;
subplot(2,1,1);plot(errorVec(:,X1),'r*');
grid on;
ylabel('X1 error');
xlabel('Sample Point');
subplot(2,1,2);plot(errorVec(:,2),'b*');
ylabel('Truck Angle error');
xlabel('Sample Point');
grid on;



