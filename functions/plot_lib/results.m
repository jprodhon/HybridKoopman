function results(EXAMPLE,eval_sets,koops,PARAMETER,do_plot)
    N = max(size(PARAMETER));
    m = size(eval_sets{1},2);
    norm_base = global_error(eval_sets{1});
    norm_fullkoop = global_error(eval_sets{1}-koops{1});
    list_normred = zeros(1,N);
    list_normredkoop = zeros(1,N);
    for i=1:N
        norm_red = global_error(eval_sets{1}-eval_sets{1+i});
        list_normred(i) = norm_red;
        norm_redkoop = global_error(...
            (eval_sets{1}-eval_sets{1+i})-koops{1+i}(:,1:m,:));
        list_normredkoop(i) = norm_redkoop;
    end
    if N == 1
        disp("dynamics norm")
        disp(norm_base)
        disp("full koopman norm")
        disp(norm_fullkoop)
        disp("reduced dynamics norm")
        disp(norm_red)
        disp("reduced koopman norm")
        disp(norm_redkoop)
    else
        figure('Visible','on')
        hold on
        plot(PARAMETER,norm_base*ones(1,N))
        plot(PARAMETER,list_normred)
        legend('full dynamics','reduced dynamics')
        xlabel('noise parameter value')
        ylabel('L2 reconstruction norm')
        title('Norm of the dynamics')
        figure('Visible','on')
        hold on
        plot(PARAMETER,norm_fullkoop*ones(1,N))
        plot(PARAMETER,list_normredkoop)
        legend('full koopman','reduced koopman')
        xlabel('noise parameter value')
        ylabel('L2 reconstruction norm')
        title('Norm of the reconstruction error')
    end
    if do_plot
        plot_sol(EXAMPLE,eval_sets,koops,PARAMETER);
    end
end

function plot_sol(EXAMPLE,eval_sets,koops,PARAMETER)
    global NX
    N = max(size(koops))-1;
    koop_recfull = koops{1};
    eval = eval_sets{1};
    for j=1:N
        koop_recred = koops{1+j}; eval_red = eval_sets{1+j};
        delete('figure_saves/*.png')
        export_sol(eval,EXAMPLE,'full-dynamics',PARAMETER(j))
        export_sol(eval-eval_red,EXAMPLE,'reduced-dynamics',PARAMETER(j))
        export_sol(koop_recfull,EXAMPLE,'koopman-on-full-dynamics',... 
            PARAMETER(j))
        export_sol(koop_recred(:,1:NX,:),EXAMPLE, ...
            'koopman-on-reduced-dynamics',PARAMETER(j))
    end
end

function export_sol(sol,EXAMPLE,type,parameter)
    global T NT L NX 
    n_dyn = size(sol,1);
    for i=1:n_dyn
        figure('Visible','on')
        surf(linspace(0,L,NX),linspace(0,T,NT),squeeze(sol(i,:,:))')
        title(sprintf('equation: %s, dynamics: %s',EXAMPLE,type))
        subtitle(sprintf('simulation number: %d, parameter: %.2f',...
            i,parameter))
        xlabel('space')
        ylabel('time')
        export_fig(sprintf('figure_saves/%.2f_%s_%d_%s.png',...
            parameter,EXAMPLE,i,type))
    end
end

function er = global_error(set)
    [n_eval,nx,nt] = size(set);
    err = zeros(1,nt);
    for n=1:nt
        error = zeros(n_eval,1);
        for m=1:n_eval
            error(m) = norm(squeeze(set(m,:,n)))/sqrt(nt*nx);
        end
        err(n) = sum(error)/n_eval;
    end
    er = sum(err);
end