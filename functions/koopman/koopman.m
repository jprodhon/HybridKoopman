function koops = koopman(eval_sets,snapshots,nSnapshots,dict)
% Generates the koopman reconstruction of the evaluation data sets
    global NT
    N = max(size(eval_sets)) - 1;
    n_eval = size(eval_sets{1},1);
    NN_init(N,dict);
    for k=1:N+1
        if k == 1
            Z = snapshots{1,1}; Zt = snapshots{1,2};
            eval = eval_sets{1};
        else
            Z = cat(1,snapshots{1,1}-snapshots{k,1},snapshots{1,1});
            Zt = cat(1,snapshots{1,2}-snapshots{k,2},snapshots{1,2});
            eval = cat(2,eval_sets{1}-eval_sets{k},eval_sets{1});
        end
        [K, error_loss] = koop_oper(Z,Zt,nSnapshots,dict,k-1);
        [eval_ini, n_func] = NN_eval_ini(k,dict);
        koop_approx = zeros(n_eval,size(Z,1),NT);
        for m=1:n_eval
            x_sol = squeeze(eval(m,:,:));
            obs_sol = observables(x_sol,dict);
            full_sol = NN_full_sol(dict,eval_ini,obs_sol,n_func,m);
            for t=2:NT; full_sol(:,t) = K*full_sol(:,t-1); end
            koop_approx(m,:,:) = full_sol(2:1+size(x_sol,1),:); 
        end
        koops{k} = koop_approx;
    end
end


function [eval_ini, n_func] = NN_eval_ini(k,dict)
% To incorporate NN dictionary in case
    if sum(contains(dict(:,1),"neural_net"))
        n_func = sum(contains(dict(:,1),"neural_net") ...
            .*cellfun(@str2num,dict(:,2)));
        if n_func > 0
            eval_ini = readmatrix(strcat('functions/dataset/saves/eval_NN',int2str(k-1),'.txt'));
        else; eval_ini = []; 
        end
    else; eval_ini = []; n_func = 0;
    end
end

function NN_init(N,dict)
% To incorporate NN dictionary in case
    if sum(contains(dict(:,1),"neural_net"))
        n_func = sum(contains(dict(:,1),"neural_net") ...
            .*cellfun(@str2num,dict(:,2)));
        if n_func > 0
            NN_dict(N,n_func,0); 
        end
    end
end


function full_sol = NN_full_sol(dict,eval_ini,obs_sol,n_func,m)
% To incorporate NN dictionary in case
    global NT
    if sum(contains(dict(:,1),"neural_net"))
        if n_func > 0
            full_sol = zeros(size(obs_sol,1)+size(eval_ini,1),NT);
            full_sol(:,1) = [obs_sol(:,1);eval_ini(:,m)];
        else; full_sol = obs_sol;
        end
    else; full_sol = obs_sol;
    end
end

function [K,error_loss] = koop_oper(Z,Zt,nSnapshots,dict,k)
% returns the Koopman operator approximation
    indices = init_for_koop(Z,nSnapshots,dict);
    Z = Z(:,indices); Zt = Zt(:,indices);
    Pz = observables(Z,dict); Pzt = observables(Zt,dict);
    if sum(contains(dict(:,1),"neural_net"))
        N = sum(contains(dict(:,1),"neural_net") ...
            .*cellfun(@str2num,dict(:,2)));
        if N > 0
            NN_obs = readmatrix(strcat('functions/dataset/saves/NN_obs',int2str(k),'.txt'));
            Pz = [Pz;NN_obs(:,indices)];
            NN_obs_t = readmatrix(strcat('functions/dataset/saves/NN_obs_t',int2str(k),'.txt'));
            Pzt = [Pzt;NN_obs_t(:,indices)];
        end
    end
    K = (Pzt*Pz')*pinv(Pz*Pz');
    error_loss = norm(Pzt-K*Pz,'fro')/sqrt(size(Pz,2));
end

function indices = init_for_koop(Z,nSnapshots,dict)
% returns the indices used for Koopman
    n_max = size(Z,2);
    indices = randperm(n_max,nSnapshots);
    N = sum(contains(dict(:,1),"RBF").*cellfun(@str2num,dict(:,2)));
    if N > 0
        [~,centers] = kmeans(Z',N,'MaxIter',500); centers = centers';
        save('functions/dataset/saves/centers.mat','centers');
    end
end

function NN_dict(N,n_func,n_neur)
% Script to run NN python script from matlab
    pyExec = 'C:\Users\Julien\anaconda3\envs\virtualenv\python.exe';
    pyRoot = fileparts(pyExec);
    p = getenv('PATH'); p = strsplit(p, ';');
    addToPath = { pyRoot
        fullfile(pyRoot, 'Library', 'mingw-w64', 'bin')
        fullfile(pyRoot, 'Library', 'usr', 'bin')
        fullfile(pyRoot, 'Library', 'bin')
        fullfile(pyRoot, 'Scripts')
        fullfile(pyRoot, 'bin') };
    p = [addToPath(:); p(:)]; p = unique(p, 'stable'); 
    p = strjoin(p, ';'); setenv('PATH', p);
    n_lay = 5;
%     param = strcat(is_full," ",int2str(nmax_sec)," ", ...
%             int2str(n_func)," ",int2str(n_lay)," ",int2str(n_neur));
    param = strcat(int2str(N)," ",int2str(n_func));
    py_path = 'functions\koopman\NN_script.py';
    command1 = 'CALL conda.bat activate virtualenv';
    command2 = strcat('python'," ",py_path," ",param);
    system(strcat(command1," && ",command2))
end