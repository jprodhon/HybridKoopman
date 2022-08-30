function nSnapshots = examples_init(EXAMPLE,N_TRAIN,nSnapshots)
% Selection of parameters for the examples
    global T NT L NX 
    if EXAMPLE=="duffing"
        T = 10; NT = 101;
    elseif EXAMPLE=="lorenz"
        T = 10; NT = 101;
    elseif EXAMPLE=="wave"
        T = 5; NT = 101; L = 1; NX = 50;
    elseif EXAMPLE=="heat"
        T = 5; NT = 101; L = 1; NX = 50; 
    elseif EXAMPLE=="KS"
        T = 10; NT = 101; L = 2*pi; NX = 50; 
    end
    nSnapshots = min(N_TRAIN*NT,nSnapshots);
end