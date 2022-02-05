%% SBLP without spike
cvx_begin
    cvx_solver gurobi_2
    variable x(num_class,num_flight,num_channel,num_period)
    variable x0(num_channel,num_period)
    expression r
    expression t0(num_flight,num_channel,num_period)
    r = price'*squeeze(sum(sum(sum(x,2),3),4));
    maximize( r )
    subject to
    x0 + squeeze(sum(sum(x,1),2)) <= arrival_rate;
    squeeze(sum(sum(sum(x,1),3),4)) <= capacity';
    for f=1:num_flight
        for j=1:num_class
            squeeze(x(j,f,:,:)./v(j,f,:,:)) <= x0./v0;
        end
    end
    x>=0; x0>=0;
cvx_end

expected_sales = squeeze(sum(sum(x,3),4));
y_sblp = cumsum(expected_sales);
z_sblp = r;

% x=full(x);
assortment_control = zeros(num_class,num_flight,num_channel,num_period);
for j=1:num_class
    for f=1:num_flight
        assortment_control(j,f,:,:) = squeeze(x(j,f,:,:)./v(j,f,:,:))./(x0./v0);
    end
end