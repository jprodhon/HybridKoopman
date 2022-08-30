function sim_sets = duffdata_set(N_SIM,parameter)
    global NT
    Tx1 = 2-4*rand(1,N_SIM); Tx2 = 2-4*rand(1,N_SIM);
    N = max(size(parameter));
    sim_sets = zeros(1+N,N_SIM,2,NT);
    for k=1:N_SIM
        x1 = Tx1(k); x2 = Tx2(k); 
        sim_sets(1,k,:,:) = solution_duff(x1,x2,1);
        for i=1:N
            sim_sets(1+i,k,:,:) = solution_duff(x1,x2,parameter(i));
        end
    end
end

function sol = solution_duff(x1,x2,parameter)
    global T NT
    lt = linspace(0,T,NT);
    a = 1; b = -1; d = parameter*0.5;
    f = @(t,x) [x(2);-d*x(2)-x(1)*(b+a*x(1)^2)];
    sol_e = ode45(f,[0 T],[x1;x2]);
    sol = deval(sol_e,lt);
end