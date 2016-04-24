function output = bipolarLog(input)
	output = (2./(1+exp(-input))) - 1;
end