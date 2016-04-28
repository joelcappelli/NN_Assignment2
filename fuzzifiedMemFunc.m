function output = fuzzifiedMemFunc(FuzzySet_Funcs,FuzzySet_output,FuzzySet_inputBreakPts,errorInput,maxMemFuncOutput)
output = {};

numFuzzySets = size(FuzzySet_output,1);
numIntersections = 0;
for i = 1:numFuzzySets
    fuzzySetoutput_temp = FuzzySet_output{i};
    fuzzySetBKpts_temp = FuzzySet_inputBreakPts{i};
    for j = 1:(size(fuzzySetoutput_temp,2)-1)
        intersection = polyval(linearEqu([fuzzySetBKpts_temp(j) fuzzySetoutput_temp(j)],[fuzzySetBKpts_temp(j+1) fuzzySetoutput_temp(j+1)]),errorInput);
        if((intersection > 0) && (intersection <= maxMemFuncOutput))
            numIntersections = numIntersections + 1;
            output(2,numIntersections) = FuzzySet_Funcs(i);
            output{1,numIntersections} = intersection;
        end
    end
end
end