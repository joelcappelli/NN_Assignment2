function output = determineRules(Fuzzified_Set1,Fuzzified_Set2,ruleType,set1_Funcs,set2_Funcs,setOutput_Funcs,FAM)

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

%rule output
numFuzzified_Set1 = size(Fuzzified_Set1,2);
numFuzzified_Set2 = size(Fuzzified_Set2,2);
output = cell(2,numFuzzified_Set1*numFuzzified_Set2);
numRules = 0;
for i = 1:numFuzzified_Set1
    for j = 1:numFuzzified_Set2
        logical_col = cellfun(cellfind(Fuzzified_Set2{2,j}),set1_Funcs);
        logical_row = cellfun(cellfind(Fuzzified_Set1{2,i}),set2_Funcs);
        findOutput = cellfun(cellfind(FAM(logical_row,logical_col)),setOutput_Funcs);
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