function [eval_sets,snaps] = data_set(EXAMPLE,N_TRAIN,N_EVAL,parameter)
% returns all the data needed for the Koopman approximation.
% the format of train_sets and eval_sets is {1} for full dynamics, {1+i} 
% for the reduced one and then for both dynamics: 
% (nb of trajectories,nx,nt)
    if EXAMPLE == "already_saved"
        [eval_sets,snaps] = load_load();
    else        
        delete('functions/dataset/saves/*')
        train_sets = get_dynamics(EXAMPLE,N_TRAIN,parameter);
        eval_sets = get_dynamics(EXAMPLE,N_EVAL,parameter);
        [snaps, eval_sets] = get_Zs(train_sets,eval_sets);
        save_save(eval_sets,snaps);
        export_ML(eval_sets,snaps);
    end
end

function sim_sets = get_dynamics(EXAMPLE,N_SIM,parameter)
    if EXAMPLE=="duffing"
        sim_sets = duffdata_set(N_SIM,parameter);
    elseif EXAMPLE=="lorenz"
        sim_sets = lorenzdata_set(N_SIM,parameter);
    elseif EXAMPLE=="wave"
        sim_sets = wavedata_set(N_SIM,parameter);
    elseif EXAMPLE=="heat"
        sim_sets = heatdata_set(N_SIM,parameter);
    elseif EXAMPLE=="KS"
        sim_sets = KSdata_set(N_SIM,parameter);
    end
end

function [snapshots, eval_sets] = get_Zs(train_sets,eval)
% Snapshots selection and eval_sets reshape
    N = size(train_sets,1)-1;
    n_eval = size(eval,2);
    [Z,Zt] = calcul_Zs(squeeze(train_sets(1,:,:,:)));
    snapshots{1,1} = Z; snapshots{1,2} = Zt; 
    if n_eval == 1
        ev = zeros(size(eval,2),size(eval,3),size(eval,4));
        ev(1,:,:) = squeeze(eval(1,1,:,:));
        eval_sets{1} = ev;
    else; eval_sets{1} = squeeze(eval(1,:,:,:)); 
    end
    for i=1:N
        [Z_red,Zt_red] = calcul_Zs(squeeze(train_sets(1+i,:,:,:)));
        snapshots{1+i,1} = Z_red; snapshots{1+i,2} = Zt_red;
        if n_eval == 1
            ev = zeros(size(eval,2),size(eval,3),size(eval,4));
            ev(1,:,:) = squeeze(eval(1+i,1,:,:));
            eval_sets{1+i} = ev;
        else; eval_sets{1+i} = squeeze(eval(1+i,:,:,:));
        end
    end
end

function [Z,Zt] = calcul_Zs(train)
    [nb_simu,m,nt] = size(train);
    Z = zeros(m,(nt-1)*nb_simu); Zt = zeros(m,(nt-1)*nb_simu);
    for i=1:nb_simu
        for j=1:m
            Z(j,(i-1)*(nt-1)+1:i*(nt-1)) = train(i,j,1:nt-1);
            Zt(j,(i-1)*(nt-1)+1:i*(nt-1)) = train(i,j,2:nt);
        end
    end
end

function [eval_sets,snapshots] = load_load()
    load('functions/dataset/saves/eval_sets.mat','eval_sets');
    load('functions/dataset/saves/snapshots.mat','snapshots');
end

function save_save(eval_sets,snapshots)
    save('functions/dataset/saves/eval_sets.mat','eval_sets');
    save('functions/dataset/saves/snapshots.mat','snapshots');
end

function export_ML(eval_sets,snapshots)
% Export data sets for NN Koopman
    N = max(size(eval_sets))-1;
    eval = eval_sets{1};
    Z = snapshots{1,1}; Zt = snapshots{1,2}; 
    writematrix(Z,'functions/dataset/saves/Z.txt')
    writematrix(Zt,'functions/dataset/saves/Zt.txt')
    eval_out = dynamics_out(eval);
    writematrix(eval_out,'functions/dataset/saves/eval.txt')
    for i=1:N
        eval_red = eval_sets{1+i};
        Z_red = snapshots{1+i,1}; Zt_red = snapshots{1+i,2};
        writematrix(Z_red,strcat('functions/dataset/saves/Z_red',int2str(i),'.txt'))
        writematrix(Zt_red,strcat('functions/dataset/saves/Zt_red',int2str(i),'.txt'))
        evalred_out = dynamics_out(cat(2,eval_red,eval));
        writematrix(evalred_out,strcat('functions/dataset/saves/eval_red',int2str(i),'.txt'))
    end
end

function dyn_out = dynamics_out(dyn)
    [n_dyn,m,NT] = size(dyn);
    dyn_out = zeros(n_dyn*m,NT);
    for i=1:n_dyn
        dyn_out((i-1)*m+1:i*m,:) = dyn(i,:,:);
    end
end