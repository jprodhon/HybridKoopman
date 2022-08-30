function sim_sets = lorenzdata_set(N_SIM,parameter)
    global NT
    N = max(size(parameter));
    Tx1 = 2*(1-2*rand(1,N_SIM)); 
    Tx2 = 2*(1-2*rand(1,N_SIM)); 
    Tx3 = 2*(1-2*rand(1,N_SIM));
    sim_sets = zeros(1+N,N_SIM,3,NT);
    for k=1:N_SIM
        x1 = Tx1(k); x2 = Tx2(k); x3 = Tx3(k);
        sim_sets(1,k,:,:) = solution_lorenz(x1,x2,x3,1);
        for i=1:N
            sim_sets(1+i,k,:,:) = solution_lorenz(x1,x2,x3,parameter(i));
        end
    end
end

function sol = solution_lorenz(x1,x2,x3,parameter)
    global T NT
    lt = linspace(0,T,NT);
    sigma = 1;
    rho = parameter*10;
    beta = 8/3;
    f = @(t,x) [sigma*(x(2)-x(1)); x(1)*(rho-x(3))-x(2); ...
        x(1)*x(2)-beta*x(3)];
    sol_e = ode45(f,[0 T],[x1;x2;x3]);
    sol = deval(sol_e,lt);
end