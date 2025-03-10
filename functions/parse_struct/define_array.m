function positions = define_array(type, baseline, ratio)

if strcmp(type, "X-Array")

    x2 = baseline / sqrt((ratio^2 + 1));
    x1 = ratio * x2;

    positions = [ x1/2,  x2/2;
                  x1/2, -x2/2;
                 -x1/2, -x2/2;
                 -x1/2,  x2/2];


elseif strcmp(type, "Linear")

    x1 = baseline / ratio;
    x2 = baseline;

    positions = [ x2/2, 0;
                  x1/2, 0;
                 -x1/2, 0;
                 -x2/2, 0];

elseif strcmp(type, "Diamond")

    x1 = baseline / ratio;
    x2 = baseline;

    positions = [ x2/2,     0;
                  0,     x1/2;
                  0,    -x1/2;
                 -x2/2,     0];

elseif strcmp(type, "Bracewell")
    
    x = baseline/ratio;

    positions = [x/2, 0,
                -x/2, 0];

else

    error("Array does not exist.")

end

end