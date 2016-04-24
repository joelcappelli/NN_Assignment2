function output = TLU(input)
    output = input;
    output(input > 0) = 1;
    output(input < 0) = -1;
end