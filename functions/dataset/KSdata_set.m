function sim_sets = KSdata_set(N_SIM,parameter)
    global NX NT L
    N = max(size(parameter));
    lx = [0:L/NX:L]'; lx = lx(1:NX);
    sim_sets = zeros(1+N,N_SIM,NX,NT);
    for k=1:N_SIM
        c1 = 0.8+rand/5; c2 = 0.5+rand/2;
        u0 = c1*sin(lx) + c2*exp(cos(lx));
        sim_sets(1,k,:,:) = dynamique_KS(u0,1);
        for i=1:N
            sim_sets(1+i,k,:,:) = dynamique_KS(u0,parameter(i));
        end
    end
end

function sol = dynamique_KS(u0,parameter)
    global NX NT T L
    nb = L/(2*pi);
    sol = zeros(NX,NT);
    sol(:,1) = u0;
    % u_t = -u*u_x - u_xx - u_xxxx, periodic boundary conditions 
    % on [0,32*pi]
    % computation is based on v = fft(u), so linear term is diagonal
    gap = 50;
    NTmax = gap*NT-1;
    v = fft(u0);
    % Precompute various scalars for ETDRK4
    h = T/(NTmax-1); % time step
    k = [0:NX/2-1 0 -NX/2+1:-1]'/(2*nb); % wave numbers
    Lf = k.^2 - k.^4; % Fourier multipliers
    E = exp(h*Lf); E2 = exp(h*Lf/2);
    M = 16; % no. of poiNTs for complex means
    r = exp(1i*pi*((1:M)-.5)/M); % roots of unity
    LR = h*Lf(:,ones(M,1)) + r(ones(NX,1),:);
    Q = h*real(mean( (exp(LR/2)-1)./LR ,2));
    f1 = h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
    f2 = h*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
    f3 = h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
    % Main loop
    g = -0.5i*k*parameter;
    for n = 1:NTmax
        Nv = g.*fft(real(ifft(v)).^2);
        a = E2.*v + Q.*Nv;
        Na = g.*fft(real(ifft(a)).^2);
        b = E2.*v + Q.*Na;
        Nb = g.*fft(real(ifft(b)).^2);
        c = E2.*a + Q.*(2*Nb-Nv);
        Nc = g.*fft(real(ifft(c)).^2);
        v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
        if mod(n,gap)==0
            u = real(ifft(v)); sol(:,n/gap+1) = u;
        end
    end
end