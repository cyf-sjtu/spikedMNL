cvx_begin
    cvx_solver gurobi_2
    variable x(num_class,num_class,num_flight,num_channel,num_period)
    variable x0(num_channel,num_period)
    expression r
    expression t0(num_flight,num_channel,num_period)
    r = squeeze(sum(sum(sum(sum(x,1),3),4),5))*price;
    maximize (r)
    x0 + squeeze(sum(sum(sum(x,1),2),3)) == arrival_rate;
    squeeze(sum(sum(sum(sum(x,1),2),4),5)) <= capacity;
    t0 = 0;
    for j = 1:num_class
        x(j,j+1:end,:,:,:) == 0;
        for k = 1:j-1
            squeeze(x(j,k,:,:,:))./squeeze(v(k,:,:,:)) <= squeeze(x(j,j,:,:,:))./squeeze(w(j,:,:,:));
        end
        t0 = t0 +  squeeze(x(j,j,:,:,:))./squeeze(w(j,:,:,:));
    end
    for f = 1:num_flight
        squeeze(t0(f,:,:)) <= x0./v0;
    end
    x>=0; x0>=0;
cvx_end

expected_sales = squeeze(sum(sum(sum(x,1),4),5));
y_sblp = cumsum(expected_sales);
assortment_control_p = zeros(num_class,num_flight,num_channel,num_period);
for j=1:num_class
    for f=1:num_flight
        assortment_control_p(j,f,:,:) = (squeeze(x(j,j,f,:,:))./squeeze(w(j,f,:,:))+squeeze(sum(x(j+1:end,j,f,:,:),1))./squeeze(v(j,f,:,:)))./(x0./v0+eps);
    end
end
assortment_control = cell(num_class,num_flight,num_channel,num_period);
for t = 1:num_period
    for c = 1:num_channel
        scale = v0(c,t)/x0(c,t);
        for f=1:num_flight
            tau = 0;
            for j=1:num_class
                for k=1:j-1
                    temp = x(j,k,f,c,t)/v(k,f,c,t)*scale;
                    if (temp>1e-6)
                        assortment_control{k,f,c,t} = [assortment_control{k,f,c,t}(:,:); tau, tau + temp];
                    end
                end
                temp = x(j,j,f,c,t)/w(j,f,c,t)*scale;
                if (temp>1e-6)
                    assortment_control{j,f,c,t} = [assortment_control{j,f,c,t}(:,:); tau, tau + temp];
                end
                tau = tau + temp;
            end
        end
    end
end
z_sblp = r;
