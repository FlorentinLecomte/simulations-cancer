function u0 = modelic(x)
    u0 = zeros(3,1);
    epsilon=0.01;
    if x <= 0.25
        u0(1) = exp(-x^2 /epsilon);
    end
    u0(2)=1-0.5*u0(1);
    u0(3)=0.5*u0(1);
end