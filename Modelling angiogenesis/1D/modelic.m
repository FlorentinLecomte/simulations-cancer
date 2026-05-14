function u0 = modelic(x)
    u0 = zeros(4,1);
    epsilon=0.01;
    if x <= 0.1
        u0(1) = exp(-x^2 /epsilon);
    end
    u0(2)=1;
    u0(3)=0;
    u0(4)=1;
end