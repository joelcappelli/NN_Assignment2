function output = checkInputRange(min,max,input)

if(input <= min)
    output = min;
elseif(input >= max)
    output = max;
else
    output = input;
end
end
