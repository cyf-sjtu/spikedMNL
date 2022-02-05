function log_lik = gam_log_lik(x,n,tau,purchase_hist)
    if(length(x)~=2*n)
        disp('Error: incompatible dimension!')
        return
    end
    v = exp(x(1:n));
    w = exp(x(n+1:end));
    log_lik = 0;
    w_sum = sum(w);
    for i = 1:tau
        choice_set = purchase_hist{i,1};
        choice = purchase_hist{i,2};
        if (choice==0)
            log_lik = log_lik ...
                + log(w_sum-sum(w(choice_set))+1) ...
                - log((sum(v(choice_set)-w(choice_set))+w_sum+1));
        else
            log_lik = log_lik ...
                + log(v(choice)) ...
                - log(sum(v(choice_set)-w(choice_set))+w_sum+1);
        end
    end
end