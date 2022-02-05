% use condor to parallelize the computation
function testing_compare_condor(idx)
load data_testing_compare;
rng(20171019+idx);

revenue_mnl = zeros(1,2);
revenue_smnl = zeros(1,2);

h=1;
% restore flight status
for f=1:num_flight
    reset(flight_mnl_standard(f));
    reset(flight_mnl_theft(f));
    reset(flight_smnl_standard(f));
    reset(flight_smnl_theft(f));
end
num_product = num_flight*num_class;
% sequence simulation begins
for t = 1:num_period
    seqArrival = generate_multichannel_arrival(num_channel,arrival_rate(:,t),1.0);
    v_ct = v(:,:,:,t);
    v_ct = reshape(v_ct,num_class*num_flight,num_channel);
    w_ct = w(:,:,:,t);
    w_ct = reshape(w_ct,num_class*num_flight,num_channel);
    for n=1:size(seqArrival,1)
        c = seqArrival(n,2);
        % mnl_standard
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_mnl_standard,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_mnl_standard(f),j);
            revenue_mnl(h,1) = revenue_mnl(h,1) + price(j);
        end
        % mnl_theft
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_mnl_theft,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_mnl_theft(f),j);
            revenue_mnl(h,2) = revenue_mnl(h,2) + price(j);
        end
        % smnl_standard
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_smnl_standard,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_smnl_standard(f),j);
            revenue_smnl(h,1) = revenue_smnl(h,1) + price(j);
        end
        % smnl_theft
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_smnl_theft,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_smnl_theft(f),j);
            revenue_smnl(h,2) = revenue_smnl(h,2) + price(j);
        end
    end
end
% output result
output=[revenue_mnl,revenue_smnl];
filename=sprintf('output/compare_h=%d.txt',idx);
dlmwrite(filename,output,'delimiter','\t','newline','pc','precision','%.1f');
end