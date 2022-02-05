num_loop = 1e2;
vec_seed = randi(1e8,1,num_loop);

for iter = 1:num_loop
seed = vec_seed(iter);
    
r_H=1; r_L=0.1; v_H=0.0; v_L=1; w_H=1; v_0=1;

tilde_v_H = 0.01;
tilde_v_L = 1;

hat_v_H = 0.01;
hat_v_L = 1;
hat_w_H = 1;

arrival_rate = 1000;
capacity = 800;
learning_rate = 0.1;

tau_1_nothing = 1;
tau_2_nothing = 1.0;
tau_1_mnl = 1;
tau_2_mnl = 1.0;
tau_1_smnl = 1;
tau_2_smnl = 1.0;
b_1_nothing = capacity;
b_2_nothing = arrival_rate*v_L/(v_H+v_L+v_0)*tau_1_nothing;
b_1_mnl = capacity;
b_2_mnl = arrival_rate*v_L/(v_H+v_L+v_0)*tau_1_mnl;
b_1_smnl = capacity;
b_2_smnl = arrival_rate*v_L/(v_H+v_L+v_0)*tau_1_smnl;

num_iter = 100;
% 1-tau_mnl,2-tau_smnl,3-b_mnl,4-b_mnl,5-tau_nothing,6-tau_nothing
vec_revenue = zeros(6,num_iter);
vec_n_H_1 = zeros(6,num_iter);
vec_n_H_2 = zeros(6,num_iter);
vec_n_L = zeros(6,num_iter);
% vec_n_0 = zeros(4,num_iter);
vec_N_HL = zeros(6,num_iter);
vec_N_H = zeros(6,num_iter);

for i=1:num_iter
    %% simulation
    seqArrival = generate_multichannel_arrival(1,arrival_rate,1.0);
    %% tau_nothing
    n_H_1 = 0;
    n_H_2 = 0;
    n_L = 0;
    n_0 = 0;
    N_HL = 0;
    N_H = 0;
    revenue = 0;
    rng(seed+i)
    for n=1:size(seqArrival,1)
        t = seqArrival(n,1);
        if (tau_1_nothing >= t)
            index = choose_product(2,v_0,[v_H;v_L],[1;1]);
            N_HL = N_HL+1;
            switch(index)
                case 1
                    n_H_1 = n_H_1 + 1;
                    revenue = revenue + r_H;
                case 2
                    n_L = n_L + 1;
                    revenue = revenue + r_L;
                case -1
                    n_0 = n_0 + 1;
                otherwise
            end
            if(n_H_1+n_H_2+n_L>=capacity), break; end
        elseif (tau_2_nothing >= t)
            index = choose_product(2,v_0,[w_H;0.0],[1;0]);
            N_H = N_H+1;
            switch(index)
                case 1
                    n_H_2 = n_H_2 + 1;
                    revenue = revenue + r_H;
                case -1
                    n_0 = n_0 + 1;
                otherwise
            end
            if(n_H_1+n_H_2+n_L>=capacity), break; end
        else
            index = 0;
        end
    end
    vec_revenue(5,i) = revenue;
    vec_n_H_1(5,i) = n_H_1;
    vec_n_H_2(5,i) = n_H_2;
    vec_n_L(5,i) = n_L;
    % vec_n_0(1,i);
    vec_N_HL(5,i) = N_HL;
    vec_N_H(5,i) = N_H;
    %% tau_mnl
    n_H_1 = 0;
    n_H_2 = 0;
    n_L = 0;
    n_0 = 0;
    N_HL = 0;
    N_H = 0;
    revenue = 0;
    rng(seed+i)
    for n=1:size(seqArrival,1)
        t = seqArrival(n,1);
        if (tau_1_mnl >= t)
            index = choose_product(2,v_0,[v_H;v_L],[1;1]);
            N_HL = N_HL+1;
            switch(index)
                case 1
                    n_H_1 = n_H_1 + 1;
                    revenue = revenue + r_H;
                case 2
                    n_L = n_L + 1;
                    revenue = revenue + r_L;
                case -1
                    n_0 = n_0 + 1;
                otherwise
            end
            if(n_H_1+n_H_2+n_L>=capacity), break; end
        elseif (tau_2_mnl >= t)
            index = choose_product(2,v_0,[w_H;0.0],[1;0]);
            N_H = N_H+1;
            switch(index)
                case 1
                    n_H_2 = n_H_2 + 1;
                    revenue = revenue + r_H;
                case -1
                    n_0 = n_0 + 1;
                otherwise
            end
            if(n_H_1+n_H_2+n_L>=capacity), break; end
        else
            index = 0;
        end
    end
    vec_revenue(1,i) = revenue;
    vec_n_H_1(1,i) = n_H_1;
    vec_n_H_2(1,i) = n_H_2;
    vec_n_L(1,i) = n_L;
    % vec_n_0(1,i);
    vec_N_HL(1,i) = N_HL;
    vec_N_H(1,i) = N_H;
    %% tau_smnl
    n_H_1 = 0;
    n_H_2 = 0;
    n_L = 0;
    n_0 = 0;
    N_HL = 0;
    N_H = 0;
    revenue = 0;
    rng(seed+i)
    for n=1:size(seqArrival,1)
        t = seqArrival(n,1);
        if (tau_1_smnl >= t)
            index = choose_product(2,v_0,[v_H;v_L],[1;1]);
            N_HL = N_HL+1;
            switch(index)
                case 1
                    n_H_1 = n_H_1 + 1;
                    revenue = revenue + r_H;
                case 2
                    n_L = n_L + 1;
                    revenue = revenue + r_L;
                case -1
                    n_0 = n_0 + 1;
                otherwise
            end
            if(n_H_1+n_H_2+n_L>=capacity), break; end
        elseif (tau_2_smnl >= t)
            index = choose_product(2,v_0,[w_H;0.0],[1;0]);
            N_H = N_H+1;
            switch(index)
                case 1
                    n_H_2 = n_H_2 + 1;
                    revenue = revenue + r_H;
                case -1
                    n_0 = n_0 + 1;
                otherwise
            end
            if(n_H_1+n_H_2+n_L>=capacity), break; end
        else
            index = 0;
        end
    end
    vec_revenue(2,i) = revenue;
    vec_n_H_1(2,i) = n_H_1;
    vec_n_H_2(2,i) = n_H_2;
    vec_n_L(2,i) = n_L;
    % vec_n_0(2,i);
    vec_N_HL(2,i) = N_HL;
    vec_N_H(2,i) = N_H;
    
%     %% b_nothing
%     n_H_1 = 0;
%     n_H_2 = 0;
%     n_L = 0;
%     n_0 = 0;
%     N_HL = 0;
%     N_H = 0;
%     revenue = 0;
%     rng(seed+i)
%     for n=1:size(seqArrival,1)
%         t = seqArrival(n,1);
%         if (n_L <= b_2_nothing && n_H_1+n_H_2+n_L<capacity)
%             index = choose_product(2,v_0,[v_H;v_L],[1;1]);
%             N_HL = N_HL+1;
%             switch(index)
%                 case 1
%                     n_H_1 = n_H_1 + 1;
%                     revenue = revenue + r_H;
%                 case 2
%                     n_L = n_L + 1;
%                     revenue = revenue + r_L;
%                 case -1
%                     n_0 = n_0 + 1;
%                 otherwise
%             end
%         elseif (n_H_1+n_H_2+n_L<capacity)
%             index = choose_product(2,v_0,[w_H;0.0],[1;0]);
%             N_H = N_H+1;
%             switch(index)
%                 case 1
%                     n_H_2 = n_H_2 + 1;
%                     revenue = revenue + r_H;
%                 case -1
%                     n_0 = n_0 + 1;
%                 otherwise
%             end
%         else
%             index = 0;
%         end
%     end
%     vec_revenue(6,i) = revenue;
%     vec_n_H_1(6,i) = n_H_1;
%     vec_n_H_2(6,i) = n_H_2;
%     vec_n_L(6,i) = n_L;
%     % vec_n_0(1,i);
%     vec_N_HL(6,i) = N_HL;
%     vec_N_H(6,i) = N_H;
%     
%     %% b_mnl
%     n_H_1 = 0;
%     n_H_2 = 0;
%     n_L = 0;
%     n_0 = 0;
%     N_HL = 0;
%     N_H = 0;
%     revenue = 0;
%     rng(seed+i)
%     for n=1:size(seqArrival,1)
%         t = seqArrival(n,1);
%         if (n_L <= b_2_mnl && n_H_1+n_H_2+n_L<capacity)
%             index = choose_product(2,v_0,[v_H;v_L],[1;1]);
%             N_HL = N_HL+1;
%             switch(index)
%                 case 1
%                     n_H_1 = n_H_1 + 1;
%                     revenue = revenue + r_H;
%                 case 2
%                     n_L = n_L + 1;
%                     revenue = revenue + r_L;
%                 case -1
%                     n_0 = n_0 + 1;
%                 otherwise
%             end
%         elseif (n_H_1+n_H_2+n_L<capacity)
%             index = choose_product(2,v_0,[w_H;0.0],[1;0]);
%             N_H = N_H+1;
%             switch(index)
%                 case 1
%                     n_H_2 = n_H_2 + 1;
%                     revenue = revenue + r_H;
%                 case -1
%                     n_0 = n_0 + 1;
%                 otherwise
%             end
%         else
%             index = 0;
%         end
%     end
%     vec_revenue(3,i) = revenue;
%     vec_n_H_1(3,i) = n_H_1;
%     vec_n_H_2(3,i) = n_H_2;
%     vec_n_L(3,i) = n_L;
%     % vec_n_0(1,i);
%     vec_N_HL(3,i) = N_HL;
%     vec_N_H(3,i) = N_H;
%     %% b_smnl
%     n_H_1 = 0;
%     n_H_2 = 0;
%     n_L = 0;
%     n_0 = 0;
%     N_HL = 0;
%     N_H = 0;
%     revenue = 0;
%     rng(seed+i)
%     for n=1:size(seqArrival,1)
%         t = seqArrival(n,1);
%         if (n_L <= b_2_smnl && n_H_1+n_H_2+n_L<capacity)
%             index = choose_product(2,v_0,[v_H;v_L],[1;1]);
%             N_HL = N_HL+1;
%             switch(index)
%                 case 1
%                     n_H_1 = n_H_1 + 1;
%                     revenue = revenue + r_H;
%                 case 2
%                     n_L = n_L + 1;
%                     revenue = revenue + r_L;
%                 case -1
%                     n_0 = n_0 + 1;
%                 otherwise
%             end
%         elseif (n_H_1+n_H_2+n_L<capacity)
%             index = choose_product(2,v_0,[w_H;0.0],[1;0]);
%             N_H = N_H+1;
%             switch(index)
%                 case 1
%                     n_H_2 = n_H_2 + 1;
%                     revenue = revenue + r_H;
%                 case -1
%                     n_0 = n_0 + 1;
%                 otherwise
%             end
%         else
%             index = 0;
%         end
%     end
%     vec_revenue(4,i) = revenue;
%     vec_n_H_1(4,i) = n_H_1;
%     vec_n_H_2(4,i) = n_H_2;
%     vec_n_L(4,i) = n_L;
%     % vec_n_0(4,i);
%     vec_N_HL(4,i) = N_HL;
%     vec_N_H(4,i) = N_H;
    %% estimate and optimize
    
    %% tau_mnl
    n_H_1 = sum(vec_n_H_1(1,:));
    n_H_2 = sum(vec_n_H_2(1,:));
    n_L = sum(vec_n_L(1,:));
    N_HL = sum(vec_N_HL(1,:));
    N_H = sum(vec_N_H(1,:));
    x = fminunc(@(x) N_HL*log(x(1)+x(2)+v_0)+N_H*log(x(1)+v_0)-(n_H_1+n_H_2)*log(x(1))-n_L*log(x(2)), [1;1]);
    tilde_v_H = x(1);
    tilde_v_L = x(2); 
    
    cvx_begin
    cvx_solver gurobi_2
    variable t_HL
    variable t_H 
    r = (r_H*tilde_v_H+r_L*tilde_v_L)/(tilde_v_H+tilde_v_L+v_0)*t_HL + ...
        (r_H*tilde_v_H)/(tilde_v_H+v_0)*t_H;
    maximize (r)
    subject to 
        t_HL + t_H <= 1;
        (tilde_v_H+tilde_v_L)/(tilde_v_H+tilde_v_L+v_0)*t_HL +...
            tilde_v_H/(tilde_v_H+v_0)*t_H <= capacity/arrival_rate;
        t_HL>=0; t_H>=0; 
    cvx_end
    tau_1_mnl = (1-learning_rate)*tau_1_mnl + learning_rate*t_HL;

    %% tau_smnl
    n_H_1 = sum(vec_n_H_1(2,:));
    n_H_2 = sum(vec_n_H_2(2,:));
    n_L = sum(vec_n_L(2,:));
    N_HL = sum(vec_N_HL(2,:));
    N_H = sum(vec_N_H(2,:));
    x = fminunc(@(x) N_HL*log(x(1)+x(2)+v_0)+N_H*log(x(3)+v_0)-n_H_1*log(x(1))-n_L*log(x(2))-n_H_2*log(x(3)), [1;1;1]);
    hat_v_H = x(1);
    hat_v_L = x(2);
    hat_w_H = x(3);
    cvx_begin
    cvx_solver gurobi_2
    variable t_HL
    variable t_H 
    r = (r_H*hat_v_H+r_L*hat_v_L)/(hat_v_H+hat_v_L+v_0)*t_HL + ...
        (r_H*hat_w_H)/(hat_w_H+v_0)*t_H;
    maximize (r)
    subject to 
        t_HL + t_H <= 1;
        (hat_v_H+hat_v_L)/(hat_v_H+hat_v_L+v_0)*t_HL +...
            hat_w_H/(hat_w_H+v_0)*t_H <= capacity/arrival_rate;
        t_HL>=0; t_H>=0; 
    cvx_end
    tau_1_smnl = (1-learning_rate)*tau_1_smnl + learning_rate*t_HL;
%     %% b_mnl
%     n_H_1 = sum(vec_n_H_1(3,:));
%     n_H_2 = sum(vec_n_H_2(3,:));
%     n_L = sum(vec_n_L(3,:));
%     N_HL = sum(vec_N_HL(3,:));
%     N_H = sum(vec_N_H(3,:));
%     x = fminunc(@(x) N_HL*log(x(1)+x(2)+v_0)+N_H*log(x(1)+v_0)-(n_H_1+n_H_2)*log(x(1))-n_L*log(x(2)), [1;1]);
%     tilde_v_H = x(1);
%     tilde_v_L = x(2);     
%     cvx_begin
%     cvx_solver gurobi_2
%     variable t_HL
%     variable t_H 
%     r = (r_H*tilde_v_H+r_L*tilde_v_L)/(tilde_v_H+tilde_v_L+v_0)*t_HL + ...
%         (r_H*tilde_v_H)/(tilde_v_H+v_0)*t_H;
%     maximize (r)
%     subject to 
%         t_HL + t_H <= 1;
%         (tilde_v_H+tilde_v_L)/(tilde_v_H+tilde_v_L+v_0)*t_HL +...
%             tilde_v_H/(tilde_v_H+v_0)*t_H <= capacity/arrival_rate;
%         t_HL>=0; t_H>=0; 
%     cvx_end
%     b_2_mnl = (1-learning_rate)*b_2_mnl + ...
%         learning_rate*t_HL*arrival_rate*tilde_v_L/(tilde_v_H+tilde_v_L+v_0)*t_HL;
%     %% b_smnl
%     n_H_1 = sum(vec_n_H_1(4,:));
%     n_H_2 = sum(vec_n_H_2(4,:));
%     n_L = sum(vec_n_L(4,:));
%     N_HL = sum(vec_N_HL(4,:));
%     N_H = sum(vec_N_H(4,:));
%     x = fminunc(@(x) N_HL*log(x(1)+x(2)+v_0)+N_H*log(x(3)+v_0)-n_H_1*log(x(1))-n_L*log(x(2))-n_H_2*log(x(3)), [1;1;1]);
%     hat_v_H = x(1);
%     hat_v_L = x(2);
%     hat_w_H = x(3);
%     cvx_begin
%     cvx_solver gurobi_2
%     variable t_HL
%     variable t_H 
%     r = (r_H*hat_v_H+r_L*hat_v_L)/(hat_v_H+hat_v_L+v_0)*t_HL + ...
%         (r_H*hat_w_H)/(hat_w_H+v_0)*t_H;
%     maximize (r)
%     subject to 
%         t_HL + t_H <= 1;
%         (hat_v_H+hat_v_L)/(hat_v_H+hat_v_L+v_0)*t_HL +...
%             hat_w_H/(hat_w_H+v_0)*t_H <= capacity/arrival_rate;
%         t_HL>=0; t_H>=0; 
%     cvx_end
%     b_2_smnl = (1-learning_rate)*b_2_smnl + ...
%         learning_rate*t_HL*arrival_rate*hat_v_L/(hat_v_H+hat_v_L+v_0)*t_HL;
end
%% output
filename=sprintf('revenue_K%d_%d.txt', ...
    capacity,iter);
dlmwrite(filename,vec_revenue, ...
    'delimiter','\t','newline','pc');
end
figure; hold on;
plot(1:num_iter,vec_revenue(1,:), '-bo')
plot(1:num_iter,vec_revenue(2,:), '--rv')
plot(1:num_iter,vec_revenue(5,:), ':md')
legend('with MNL', 'with spiked-MNL','do nothing')
xlabel('iteration')
ylabel('revenue')

% figure; hold on;
% plot(1:num_iter,vec_revenue(3,:), '-bo')
% plot(1:num_iter,vec_revenue(4,:), '--rv')
% plot(1:num_iter,vec_revenue(6,:), ':md')
% legend('with MNL', 'with spiked-MNL','do nothing')
% xlabel('iteration')
% ylabel('revenue')

