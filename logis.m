function output = logis(input,deriv)
    if((nargin == 2) && strcmp(deriv,'deriv'))
        output = input.*(1-input);
    else
        output = 1./(1+exp(-input));
    end
end