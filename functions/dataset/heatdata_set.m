function sim_sets = heatdata_set(N_SIM,parameter)
    global NT NX la r
    N = max(size(parameter));
    sim_sets = zeros(1+N,N_SIM,NX,NT);
    for k=1:N_SIM
        la = rand; r = rand;
        sim_sets(1,k,:,:) = dynamique_heat(1);
        for i=1:N
            sim_sets(1+i,k,:,:) = dynamique_heat(parameter(i));
        end
    end
end

function sol = dynamique_heat(parameter)
    global T NT L NX p
    p = parameter;
    m = 0;
    x = linspace(0,L,NX);
    t = linspace(0,T,NT);
    sol = pdepe(m,@heatpde,@heatic,@heatbc,x,t)';
end

function [c,f,s] = heatpde(x,t,u,dudx)
    global p
    c = p*25; f = dudx; s = 1;
end

function u0 = heatic(x)
    global r la
    lg = 0.1+0.2*la;
    deb = 0.1+(1-lg-0.2)*r;
    if x>deb && x<deb+lg
        u0 = 1;
    else
        u0 = 0;
    end
end

function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
    pl = ul; ql = 0; pr = ur; qr = 0;
end