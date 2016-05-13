function [output, modifiedFAM] = determineRules(Fuzzified_Set1,Fuzzified_Set2,ruleType,set1_Funcs,set2_Funcs,setOutput_Funcs,FAM,findDominant)

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
modifiedFAM = FAM;

%rule output
numFuzzified_Set1 = size(Fuzzified_Set1,2);
numFuzzified_Set2 = size(Fuzzified_Set2,2);
output = cell(2,numFuzzified_Set1*numFuzzified_Set2);
numRules = 0;
dominantRuleFuzzySet1 = 0;
dominantRuleFuzzySet2 = 0;
dominantRuleFireStrength = 0;

for i = 1:numFuzzified_Set1
    for j = 1:numFuzzified_Set2
        logical_col = cellfun(cellfind(Fuzzified_Set2{2,j}),set2_Funcs);
        logical_row = cellfun(cellfind(Fuzzified_Set1{2,i}),set1_Funcs);
        findOutput = cellfun(cellfind(FAM(logical_row,logical_col)),setOutput_Funcs);
        RuleOutputFunc = setOutput_Funcs(findOutput);
        numRules = numRules + 1;
        
        if(strcmp(ruleType,'MOM-prod'))
            output{1,numRules} = Fuzzified_Set1{1,i}*Fuzzified_Set2{1,j};
        elseif(strcmp(ruleType,'COA-min'))
            output{1,numRules} = min(Fuzzified_Set1{1,i},Fuzzified_Set2{1,j});
        end
        
        output{2,numRules} = RuleOutputFunc{:};
        
        if(output{1,numRules} > dominantRuleFireStrength)
            dominantRuleFireStrength = output{1,numRules};
            dominantRuleFuzzySet1 = i;
            dominantRuleFuzzySet2 = j;
        end
    end
end

weaker = 1;

if(findDominant)
    logical_col = cellfun(cellfind(Fuzzified_Set2{2,dominantRuleFuzzySet2}),set2_Funcs);
    logical_row = cellfun(cellfind(Fuzzified_Set1{2,dominantRuleFuzzySet1}),set1_Funcs);
    findOutput = cellfun(cellfind(FAM(logical_row,logical_col)),setOutput_Funcs);
    if(weaker)
        findOutput = [0 findOutput(1:end-1)];
        if(isempty(find(findOutput)))
            findOutput(end) = 1;
        end
    else
        findOutput = [findOutput(2:end) 0];
        if(isempty(find(findOutput)))
            findOutput(1) = 1;
        end
    end
    
    dominantRuleOutputFunc = setOutput_Funcs(logical(findOutput));

    modifiedFAM(logical_row,logical_col) = dominantRuleOutputFunc;
    numRules = 0;    

    for i = 1:numFuzzified_Set1
        for j = 1:numFuzzified_Set2
            logical_col = cellfun(cellfind(Fuzzified_Set2{2,j}),set2_Funcs);
            logical_row = cellfun(cellfind(Fuzzified_Set1{2,i}),set1_Funcs);
            findOutput = cellfun(cellfind(modifiedFAM(logical_row,logical_col)),setOutput_Funcs);
            RuleOutputFunc = setOutput_Funcs(findOutput);
            numRules = numRules + 1;

            if(strcmp(ruleType,'MOM-prod'))
                output{1,numRules} = Fuzzified_Set1{1,i}*Fuzzified_Set2{1,j};
            elseif(strcmp(ruleType,'COA-min'))
                output{1,numRules} = min(Fuzzified_Set1{1,i},Fuzzified_Set2{1,j});
            end

            output{2,numRules} = RuleOutputFunc{:};
        end
    end
end

end