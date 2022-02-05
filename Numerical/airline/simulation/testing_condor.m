% use condor to parallelize the computation
function testing_condor(idx)
load data_post_training;
rng(20171019+idx);

revenue_fcfs = zeros(1,1);
revenue_robust = zeros(1,2);
revenue_EMSRb = zeros(1,2);
revenue_updated = zeros(1,2);
revenue_sblp = zeros(1,3);
%revenue_kw = zeros(1,2);

h=1;
% restore flight status
for f=1:num_flight
    reset(flight_fcfs(f));
    reset(flight_robust_standard(f));
    reset(flight_robust_theft(f));
    %reset(flight_robust_dynamic_standard(f));
    %reset(flight_robust_dynamic_theft(f));
    reset(flight_EMSRb_standard(f));
    reset(flight_EMSRb_theft(f));
    reset(flight_sblp_standard(f));
    reset(flight_sblp_theft(f));
    reset(flight_sblp_assortment(f));
    reset(flight_updated_standard(f));
    reset(flight_updated_theft(f));
    %reset(flight_kw_standard(f));
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
        % current_time = time + seqArrival(n);

        % control policy specific
        % FCFS
        % spike = get_spike(num_flight,num_class,flight_fcfs);
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_fcfs,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_fcfs(f),j);
            revenue_fcfs(h,1) = revenue_fcfs(h,1) + price(j);
        end
        % robust_standard
        % spike = get_spike(num_flight,num_class,flight_robust_standard);
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_robust_standard,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_robust_standard(f),j);
            revenue_robust(h,1) = revenue_robust(h,1) + price(j);
        end
        % robust_theft
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_robust_theft,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_robust_theft(f),j);
            revenue_robust(h,2) = revenue_robust(h,2) + price(j);
        end
        % EMSRb_standard
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_EMSRb_standard,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_EMSRb_standard(f),j);
            revenue_EMSRb(h,1) = revenue_EMSRb(h,1) + price(j);
        end
        % EMSRb_theft
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_EMSRb_theft,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_EMSRb_theft(f),j);
            revenue_EMSRb(h,2) = revenue_EMSRb(h,2) + price(j);
        end
        % sblp_standard
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_sblp_standard,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_sblp_standard(f),j);
            revenue_sblp(h,1) = revenue_sblp(h,1) + price(j);
        end
        % sblp_theft
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_sblp_theft,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_sblp_theft(f),j);
            revenue_sblp(h,2) = revenue_sblp(h,2) + price(j);
        end
        % sblp_assortment
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_sblp_assortment,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_sblp_assortment(f),j);
            revenue_sblp(h,3) = revenue_sblp(h,3) + price(j);
        end
        % updated_standard
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_updated_standard,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_updated_standard(f),j);
            revenue_updated(h,1) = revenue_updated(h,1) + price(j);
        end
        % updated_theft
        [avail,spike] = get_avail_spike(num_flight,num_class,flight_updated_theft,c,t);
        vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
        index = choose_product(num_product,v0(c,t),vtemp,avail);
        % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
        if(index > 0)
            [f,j] = find_flight_product(index,num_flight,num_class);
            put_booking(flight_updated_theft(f),j);
            revenue_updated(h,2) = revenue_updated(h,2) + price(j);
        end
        % KW_standard
%         [avail,spike] = get_avail_spike(num_flight,num_class,flight_kw_standard,c,t);
%         vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
%         index = choose_product(num_product,v0(c,t),vtemp,avail);
%         % [U, I] = generate_customer(Vtemp,v0(c,t),seed);
%         if(index > 0)
%             [f,j] = find_flight_product(index,num_flight,num_class);
%             put_booking(flight_kw_standard(f),j);
%             revenue_kw(h,1) = revenue_kw(h,1) + price(j);
%         end
    end
end
% output result
output=[revenue_fcfs,...
        revenue_robust,...
        revenue_EMSRb,...
        revenue_sblp,...
        revenue_updated];
filename=sprintf('output_training/h=%d.txt',idx);
dlmwrite(filename,output,'delimiter','\t','newline','pc','precision','%.1f');
end