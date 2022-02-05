% meta-parameter
n = 15; % # of products
m = 30; % # of customer types
tau = .5e3; % training data history length
test_tau = 3e3; % testing data history length;
% drop_prob = 0.2; % when generating a customer type, drop an ite m
% shuffle_prob = 0.2; % when generating a customer type, shuffle two adjacent items
% choice_set_drop_prob = 0.2; % when generating a choice set, drop an item
eps = 1e-6;

num_ground_choice_model = 20;
num_revenue_sample = 20;
% storage for results
% columns: MNL, SMNL
vec_log_lik_train = zeros(num_ground_choice_model,3); 
vec_log_lik_test = zeros(num_ground_choice_model,3); 
vec_exitflag = zeros(num_ground_choice_model,3); 
vec_revenue = zeros(num_ground_choice_model,num_revenue_sample,4); % add offering everything
%vec_avg_rev = zeros(num_ground_choice_model,3);
vec_beats = zeros(num_ground_choice_model,4,4);
for sample_idx = 1:num_ground_choice_model
    sample_idx
% set seed
seed = 19260817 + sample_idx*6;
rng(seed)

%% generate customer types 
customer_type_list = cell(m,1);
temp_list = zeros(m,1);
iter = 1;
while(iter<=m)
    % cus_len = randi([2,5]);
    % cus_len = 5;
    % j_top = randi(round(n/2)); % highest affordable product
    j_top = randi(n);
    % j_bottom = randi([j_top,min(n,j_top+cus_len)]);  % lowest acceptable product
    j_bottom = randi([j_top,n]);
    temp = j_top:j_bottom;
    len = length(temp);
%     if (len>=2) % drop 
%         temp = temp(rand(1,len)>drop_prob);
%     end
%     len = length(temp);
%     if (len == 0)
%         continue; % skip trivial customer type
%     elseif (len>=2 && rand()>shuffle_prob) % shuffle
%         idx = randi(len-1);
%         temp([idx+1,idx]) = temp([idx,idx+1]);
%         idx = randi(len-1);
%         temp([idx+1,idx]) = temp([idx,idx+1]);
%     end
    temp = fliplr(temp);
    temp_abbr = sum(temp.*(n+1).^(0:1:len-1));
    if(sum(temp_abbr==temp_list)==0) 
        customer_type_list{iter,1} = temp;
        temp_list(iter) = temp_abbr;
        iter = iter + 1;
    else
        continue; % skip repeated customer type
    end
end

%% generate history
purchase_hist = cell(tau,2);
iter = 1;
whole_set = 1:n;
Y = zeros(tau,1);
X = zeros(tau,m);
while iter <= tau
    % choice_set = whole_set(rand(1,n)>rand());
%     while(1)
%         choice_set = 1:randi(n);
%         choice_set = choice_set(rand(1,length(choice_set))>choice_set_drop_prob);
%         if(~isempty(choice_set))
%             break
%         end
%     end
    choice_set = 1:randi(n);
    cus_type = randi(m);   
    temp = customer_type_list{cus_type,1};
    idx = find(sum(temp==choice_set',1));
    if(isempty(idx))
        choice = 0;
    else
        choice = temp(idx(1));
    end
    purchase_hist{iter,1} = choice_set;
    purchase_hist{iter,2} = choice;
    Y(iter) = choice;
    X(iter,choice_set) = 1;
    iter = iter + 1;
end

%% calibrate choice models

% MNL
options = optimoptions('fminunc','MaxFunctionEvaluations',1e4);
[x,fval,exitflag,output] = fminunc(@(x) -mnl_log_lik(x,n,tau,purchase_hist),...
                                    zeros(n,1),options);
v_mnl = exp(x);
vec_log_lik_train(sample_idx,1) = -fval;
vec_exitflag(sample_idx,1) = exitflag;
% SMNL
options = optimoptions('fminunc','MaxFunctionEvaluations',1e4);
[x,fval,exitflag,output] = fminunc(@(x) -smnl_log_lik(x,n,tau,purchase_hist),...
                                    zeros(2*n-1,1),options);
% [x,fval,exitflag,output] = fmincon(@(x) -smnl_log_lik(x,n,tau,purchase_hist),...
%                     [zeros(n,1);ones(n,1)],...
%                     [eye(n),-eye(n)],zeros(n,1),...
%                     [zeros(1,n-1),1,zeros(1,n-1),-1],0,[],[],[],options);
v_smnl = exp(x(1:n));
w_smnl = exp(x(n+1:end));
vec_log_lik_train(sample_idx,2) = -fval;
vec_exitflag(sample_idx,2) = exitflag;

% GAM
options = optimoptions('fminunc','MaxFunctionEvaluations',1e4);
[x,fval,exitflag,output] = fminunc(@(x) -gam_log_lik(x,n,tau,purchase_hist),...
                                    zeros(2*n,1),options);
v_gam = exp(x(1:n));
w_gam = exp(x(n+1:end));
vec_log_lik_train(sample_idx,3) = -fval;
vec_exitflag(sample_idx,3) = exitflag;

%% predict customer purchases (log-likelihood)
testing_data = cell(test_tau,2);
iter = 1;
whole_set = 1:n;
Xtext = zeros(test_tau,m);
Ytext = zeros(test_tau,m);
while iter <= test_tau
    %choice_set = whole_set(rand(1,n)>rand());
    choice_set = 1:randi(n);
%     if(isempty(choice_set))
%         continue
%     end
    cus_type = randi(m);
    temp = customer_type_list{cus_type,1};
    idx = find(sum(temp==choice_set',1));
    if(isempty(idx))
        choice = 0;
    else
        choice = temp(idx(1));
    end
    testing_data{iter,1} = choice_set;
    testing_data{iter,2} = choice;
    iter = iter + 1;
end
vec_log_lik_test(sample_idx,1) = mnl_log_lik(log(v_mnl),n,test_tau,testing_data);
vec_log_lik_test(sample_idx,3) = smnl_log_lik(log([v_smnl;w_smnl]),n,test_tau,testing_data);
vec_log_lik_test(sample_idx,2) = gam_log_lik(log([v_gam;w_gam]),n,test_tau,testing_data);

%% find optimal assortment
for k = 1:num_revenue_sample
    r = sort(rand(n,1)*100,'descend'); % price vector; high to low
    whole_set = 1:n;
    % OPT
    cvx_begin
        cvx_solver gurobi_2
        variable x(n,1) binary
        variable y(n,m)
        % variable y0(1,m)
        obj = sum(r'*y);
        maximize obj
        subject to
            sum(y,1) <= ones(1,m);
            y <= x*ones(1,m);
            for cus_type = 1:m
                temp = customer_type_list{cus_type,1};
                len = length(temp);
                for i = 1:len
                    for j = i+1:len
                        1 - x(temp(i),1) >= y(temp(j),cus_type); 
                    end
                end
                temp = setdiff(whole_set,temp);
                for i = temp
                    y(i,cus_type) == 0;
                end
            end
    cvx_end
    S_opt = whole_set(x>eps);
    % MNL
    cvx_begin
        cvx_solver gurobi_2
        variable x(n,1) nonnegative
        variable x0 nonnegative
        obj = r'*x;
        maximize obj
        subject to 
            sum(x)+x0 == 1;
            x./v_mnl <= x0;
    cvx_end
    S_mnl = whole_set(x>eps);
    % smnl
    a = zeros(n,1);
    b = zeros(n,1);
    for i = 1:n-1
        a(i) = 1 + sum(v_smnl(1:i-1))/w_smnl(i);
        b(i) = r(i) + sum(r(1:i-1).*v_smnl(1:i-1))/w_smnl(i);
    end
    a(n) = 1 + sum(v_smnl(1:n-1))/v_smnl(n);
    b(n) = r(i) + sum(r(1:n-1).*v_smnl(1:n-1))/v_smnl(n);
    cvx_begin
        cvx_solver gurobi_2
        variable x0 nonnegative
        variable x(n,1) nonnegative
        obj = b'*x;
        maximize obj
        subject to
            sum(a.*x)+x0 == 1;
            sum(x(1:n-1)./w_smnl)+x(n)/v_smnl(n) <= x0;
    cvx_end
    temp = find(x>eps);
    S_smnl = 1:temp;
    % gam
    cvx_begin
        cvx_solver gurobi_2
        variable x(n,1) nonnegative
        variable x0 nonnegative
        obj = r'*x;
        maximize obj
        subject to 
            sum((v_gam-w_gam)./v_gam.*x)+(1+sum(w_gam))*x0 == 1;
            x./v_gam <= x0;
    cvx_end
    S_gam = whole_set(x>eps);
    % Compute expected revenue
    revenue_opt = 0;
    revenue_mnl = 0;
    revenue_smnl = 0;
    revenue_gam = 0;
    for cus_type = 1:m % equal probability
        temp = customer_type_list{cus_type,1};
        % OPT
        idx = find(sum(temp==S_opt',1));
        if(~isempty(idx)) % if buying something
            choice = temp(idx(1));
            revenue_opt = revenue_opt + r(choice);
        end
        % MNL
        idx = find(sum(temp==S_mnl',1));
        if(~isempty(idx)) % if buying something
            choice = temp(idx(1));
            revenue_mnl = revenue_mnl + r(choice);
        end
        % SMNL
        idx = find(sum(temp==S_smnl',1));
        if(~isempty(idx)) % if buying something
            choice = temp(idx(1));
            revenue_smnl = revenue_smnl + r(choice);
        end       
        % GAM
        idx = find(sum(temp==S_gam',1));
        if(~isempty(idx)) % if buying something
            choice = temp(idx(1));
            revenue_gam = revenue_gam + r(choice);
        end       
    end
    vec_revenue(sample_idx,k,1) = revenue_mnl;
    vec_revenue(sample_idx,k,3) = revenue_smnl;
    vec_revenue(sample_idx,k,2) = revenue_gam;
    vec_revenue(sample_idx,k,4) = revenue_opt;
end


end

vec_ratio = vec_revenue(:,:,1:3)./vec_revenue(:,:,4);
vec_avg_ratio = squeeze(mean(vec_ratio,2));
vec_avg_rev = squeeze(mean(vec_revenue,2))/1000;

for iter = 1:num_ground_choice_model
for k = 1:num_revenue_sample
vec_beats(iter,:,:) = squeeze(vec_beats(iter,:,:)) +...
(squeeze(vec_revenue(iter,k,:))'>squeeze(vec_revenue(iter,k,:)));
end
end

ratio_mnl = vec_ratio(:,:,1);
ratio_mnl = ratio_mnl(:);
se_mnl=std(ratio_mnl)/sqrt(num_ground_choice_model*num_revenue_sample)
ratio_gam = vec_ratio(:,:,2);
ratio_gam = ratio_gam(:);
se_gam=std(ratio_gam)/sqrt(num_ground_choice_model*num_revenue_sample)
ratio_smnl = vec_ratio(:,:,3);
ratio_smnl = ratio_smnl(:);
se_smnl=std(ratio_smnl)/sqrt(num_ground_choice_model*num_revenue_sample)
%%
num_product = 8;
S = 1:num_product;
ct = zeros(num_product,1);
for cus_type = 1:m % equal probability
    temp = customer_type_list{cus_type,1};
    % OPT
    idx = find(sum(temp==S',1));
    if(~isempty(idx)) % if buying something
        choice = temp(idx(1));
        ct(choice)=ct(choice)+1;
    end
end
subplot(2,2,3);plot(S',flip(ct)/sum(ct));
xlim([0,num_product+1]);
ylim([0,0.8]);
xticks(1:num_product);xticklabels(num_product:-1:1)
xlabel('Product');ylabel('Fraction of Purchases');
title('n=15, p=40')