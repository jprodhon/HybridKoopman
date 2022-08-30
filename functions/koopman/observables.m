function psi = observables(data,dict)
% returns the observable vector/matrix from data with respect to the
% dictionary dict.
    [~,NT] = size(data);
    psi = [];
    if sum(contains(dict(:,1),"const").*cellfun(@str2num,dict(:,2))) > 0
        psi = [ones(1,NT); psi];
    end
    if sum(contains(dict(:,1),"state").*cellfun(@str2num,dict(:,2))) > 0
        psi = [psi; data];
    end
    if sum(contains(dict(:,1),"derivt").*cellfun(@str2num,dict(:,2))) > 0
        psi = [psi; deriv_t(data)];
    end
    if sum(contains(dict(:,1),"derivx").*cellfun(@str2num,dict(:,2))) > 0
        psi = [psi; deriv_x(data)];
    end
    if sum(contains(dict(:,1),"derivt2").*cellfun(@str2num,dict(:,2))) > 0
        psi = [psi; deriv_t(deriv_t(data))];
    end
    if sum(contains(dict(:,1),"derivx2").*cellfun(@str2num,dict(:,2))) > 0
        psi = [psi; deriv_x(deriv_x(data))];
    end
    N = sum(contains(dict(:,1),"mono").*cellfun(@str2num,dict(:,2)));
    if N > 0
        psi = [psi;psi_mono(psi,N)];
    end
    N = sum(contains(dict(:,1),"hermite").*cellfun(@str2num,dict(:,2)));
    if N > 0
        psi = [psi;psi_hermite(psi,N)];
    end
    N = sum(contains(dict(:,1),"RBF").*cellfun(@str2num,dict(:,2)));
    if N > 0
        psi = [psi;psi_RBF(data,N)];
    end
end

function psi = psi_mono(data,N)
    [m,NT] = size(data);
    if m<4
        expo = exposants(m,N);
        NN = size(expo,1);
        psi = ones(NN,NT);
        for k=1:NN
            for i=1:m
                ex = expo(k,i); 
                psi(k,:) = psi(k,:).*(data(i,:).^ex);
            end
        end
    else
        global NX
        mm = floor(m/NX);
        ind = 1; exp = 2;
        psi = zeros(N*NX,NT);
        for k=1:N
            psi((k-1)*NX+1:k*NX,:) = data((ind-1)*NX+1:ind*NX,:).^exp;
            ind = ind+1;
            if ind > mm; ind = 1; exp = exp+1; end
        end
    end
end

function psi = psi_hermite(data,N)
    [m,NT] = size(data);
    A = ones(m,N,NT); A(:,2,:) = data;
    for k=2:N-1
        A(:,k+1,:) = data.*squeeze(A(:,k,:)) - (k-1).*squeeze(A(:,k-1,:));
    end
    if m<4
        expo = exposants(m,N);
        psi = ones(N,NT); psi(2:1+m,:) = data;
        for k=2+m:N
            for i=1:m
                ex = expo(k,i)+1;
                psi(k,:) = psi(k,:).*squeeze(A(i,ex,:))';
            end
        end
    else 
        psi = ones(N*m,NT);
        for k=1:N
            psi((k-1)*m+1:k*m,:) = squeeze(A(:,k,:));
        end
    end
end

function psi = psi_RBF(data,N)
    load('functions/dataset/saves/centers.mat','centers');
    [~,NT] = size(data);
    psi = ones(N,NT);
    for k=1:N
        r = vecnorm(data-centers(:,k));
        psi(k,:) = r.^2.*log(r+0.0001);
    end
end

function expo = exposants(m,N)
    L = [zeros(1,m)];
    for i=1:N
        L = [L,i*ones(1,m)];
    end
    a = unique(nchoosek(L,m),'rows');
    nn = size(a,1);
    a = sort_tab(a,m);
    expo = [];
    for k=1:nn
        expo = [expo;unique(perms(a(k,:)),'rows')];
    end
    expo = sort_tab(expo,m);
end

function x = sort_tab(x,m)
    s = size(x,1);
    for_sort = zeros(s,m+1); for_sort(:,2:m+1) = x;
    for i=1:m
        for_sort(:,1) = for_sort(:,1) + (i^1.01).*x(:,i).^2;
    end
    sorted = sortrows(for_sort); x = sorted(:,2:m+1);
end

function der_t = deriv_t(A)
    global T
    [m,NT] = size(A);
    dt = T/(NT-1);
    der_t = zeros(m,NT);
    for i=1:3
        der_t(:,i) = -49*A(:,i)/20 + 6*A(:,i+1) - 15*A(:,i+2)/2 ...
            + 20*A(:,i+3)/3 - 15*A(:,i+4)/4 + 6*A(:,i+5)/5 - A(:,i+6)/6;
    end
    for i=4:NT-3
        der_t(:,i) = -A(:,i-3)/60 + 3*A(:,i-2)/20 - 3*A(:,i-1)/4 ...
            + 3*A(:,i+1)/4 - 3*A(:,i+2)/20 + A(:,i+3)/60;
    end
    for i=NT-2:NT
        der_t(:,i) = 49*A(:,i)/20 - 6*A(:,i-1) + 15*A(:,i-2)/2 ...
            - 20*A(:,i-3)/3 + 15*A(:,i-4)/4 - 6*A(:,i-5)/5 + A(:,i-6)/6;
    end
    der_t = der_t/dt;
end

function der_x = deriv_x(A)
    global L NX
    dx = L/(NX-1);
    [m,NT] = size(A);
    der_x = zeros(m,NT);
    for j=1:3
        der_x(j,:) = -49*A(j,:)/20 + 6*A(j+1,:) - 15*A(j+2,:)/2 ...
            + 20*A(j+3,:)/3 - 15*A(j+4,:)/4 + 6*A(j+5,:)/5 - A(j+6,:)/6;
    end
    for j=4:m-3
        der_x(j,:) = -A(j-3,:)/60 + 3*A(j-2,:)/20 - 3*A(j-1,:)/4 ...
            + 3*A(j+1,:)/4 - 3*A(j+2,:)/20 + A(j+3,:)/60;
    end
    for j=m-2:m
        der_x(j,:) = 49*A(j,:)/20 - 6*A(j-1,:) + 15*A(j-2,:)/2 ...
            - 20*A(j-3,:)/3 + 15*A(j-4,:)/4 - 6*A(j-5,:)/5 + A(j-6,:)/6;
    end
    der_x = der_x/dx;
end