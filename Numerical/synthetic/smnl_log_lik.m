function log_lik = smnl_log_lik(x,n,tau,purchase_hist)
    if(length(x)~=2*n-1)
        disp('Error: incompatible dimension!')
        return
    end
    v = exp(x(1:n));
    w = exp(x(n+1:end));
    log_lik = 0;
    for i = 1:tau
        choice_set = purchase_hist{i,1};
        choice = purchase_hist{i,2};
        top_choice = choice_set(end);
        if (choice==0)
            if (sum(choice_set==n)==0)
                log_lik = log_lik - log(sum(v(choice_set))-v(top_choice)+w(top_choice)+1);
            else
                log_lik = log_lik - log(sum(v(choice_set))+1);
            end
        elseif (choice==n)
            log_lik = log_lik + log(v(choice)) - log(sum(v(choice_set))+1);
        elseif (choice==top_choice)
            log_lik = log_lik + log(w(choice)) - log(sum(v(choice_set))-v(top_choice)+w(top_choice)+1);
        else
            if (sum(choice_set==n)==0)
                log_lik = log_lik + log(v(choice)) - log(sum(v(choice_set))-v(top_choice)+w(top_choice)+1);
            else
                log_lik = log_lik + log(v(choice)) - log(sum(v(choice_set))+1);
            end
        end
    end
end