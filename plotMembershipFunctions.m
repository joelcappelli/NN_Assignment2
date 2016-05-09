function plotMembershipFunctions(contOutputFuzzySet_output,contOutputFuzzySet_inputBreakPts,setOutput,contOutputTextXpos,...
                                xPosErrorFuzzySet_output, xPosErrorFuzzySet_inputBreakPts,setXPosError,errorTextXpos,...
                                TruckAngleErrorFuzzySet_output, TruckAngleErrorFuzzySet_inputBreakPts,setTruckAngleError,truckAngleErrorTextXpos,...
                                initialXPosError,initialTruckAngleError)

numxPosErrorFuzzySets = size(xPosErrorFuzzySet_output,1);
numTruckAngleErrorFuzzySets = size(TruckAngleErrorFuzzySet_output,1);                           
numContOutputFuzzySets = size(contOutputFuzzySet_output,1);

figure;
for i = 1:numTruckAngleErrorFuzzySets
str = setTruckAngleError{i};
subplot(3,1,1); plot(TruckAngleErrorFuzzySet_inputBreakPts{i},TruckAngleErrorFuzzySet_output{i},'b');
subplot(3,1,1); text(truckAngleErrorTextXpos(i),1,str);
hold on;
end
subplot(3,1,1);plot([initialTruckAngleError initialTruckAngleError],[0 1],'--r','linewidth',2);
grid on;
ylabel('Truck Angle (\phi deg)');

for i = 1:numxPosErrorFuzzySets
str = setXPosError{i};
subplot(3,1,2); plot(xPosErrorFuzzySet_inputBreakPts{i},xPosErrorFuzzySet_output{i},'b');
subplot(3,1,2);text(errorTextXpos(i),1,str);
hold on;
end 
subplot(3,1,2);plot([initialXPosError initialXPosError],[0 1],'--k','linewidth',2);
grid on;
ylabel('Position (x metres)');

for i = 1:numContOutputFuzzySets
str = setOutput{i};
subplot(3,1,3); plot(contOutputFuzzySet_inputBreakPts{i},contOutputFuzzySet_output{i},'b');
subplot(3,1,3);text(contOutputTextXpos(i),1,str);
hold on;
end
grid on;
ylabel('Control Angle (\theta deg)');
 
                            
end