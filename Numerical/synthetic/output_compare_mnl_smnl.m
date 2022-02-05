num_loop = 200;
num_iter = 100;
revenue_tau_nothing = zeros(num_loop,num_iter);
revenue_tau_mnl = zeros(num_loop,num_iter);
revenue_tau_smnl = zeros(num_loop,num_iter);
revenue_b_nothing = zeros(num_loop,num_iter);
revenue_b_mnl = zeros(num_loop,num_iter);
revenue_b_smnl = zeros(num_loop,num_iter);

capacity = 800;

for loop = 1:num_loop
    filename=sprintf('revenue_K%d_%d.txt', capacity,loop);
    fileID = fopen(filename);
    C = textscan(fileID,repmat('%d ',1,num_iter),'Delimiter','\t');
    fclose(fileID);
    for iter = 1:num_iter
        revenue_tau_mnl(loop,iter) = C{1,iter}(1);
        revenue_tau_smnl(loop,iter) = C{1,iter}(2);
        revenue_b_mnl(loop,iter) = C{1,iter}(3);
        revenue_b_smnl(loop,iter) = C{1,iter}(4);
        revenue_tau_nothing(loop,iter) = C{1,iter}(5);
        revenue_b_nothing(loop,iter) = C{1,iter}(6);
    end
end

% figure; hold on;
% subplot(1,2,2);hold on;
% plot(1:num_iter,mean(revenue_tau_mnl), '-bo')
% plot(1:num_iter,mean(revenue_tau_smnl), '--rv')
% plot(1:num_iter,2000*ones(1,num_iter), '-m')
% legend('with MNL', 'with spiked-MNL','z^{CDLP}')
% xlabel('Iteration')
% ylabel('Average Revenue')
% title(sprintf('Capacity=%d', capacity))

figure; hold on;
plot(1:num_iter,mean(revenue_b_mnl), '-bo')
plot(1:num_iter,mean(revenue_b_smnl), '--rv')
% plot(1:num_iter,mean(revenue_b_nothing), ':md')
plot(1:num_iter,2000*ones(1,num_iter), '-m')
legend('MNL', 'Spiked-MNL','Upper Bound')
xlabel('Iteration','fontsize',16)
ylabel('Average Revenue','fontsize',16)
title('Revenue Evolution','fontsize',20)
% title(sprintf('Booking Limit, Capacity=%d', capacity))