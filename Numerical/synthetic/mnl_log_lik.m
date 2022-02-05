function log_lik = mnl_log_lik(U,n,tau,purchase_hist)
    if(length(U)~=n)
        disp('Error: incompatible dimension!')
        return
    end
    v = exp(U);
    log_lik = 0;
    for i = 1:tau
        choice_set = purchase_hist{i,1};
        choice = purchase_hist{i,2};
        if (choice==0)
            log_lik = log_lik - log((sum(v(choice_set))+1));
        else
            log_lik = log_lik + log(v(choice)) - log(sum(v(choice_set))+1);
        end
    end
end