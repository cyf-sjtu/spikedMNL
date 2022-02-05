% FCFS
for f=1:num_flight
    flight_fcfs(f) = Flight(num_class, price, zeros(num_class-1,1),capacity(f),'standard','booking',[]);
end

% Robust
for f=1:num_flight
    y = get_robust_protection_level(num_class, price,capacity(f));
    flight_robust_standard(f) = Flight(num_class, price, y, capacity(f),'standard','booking',[]);
    flight_robust_theft(f) = Flight(num_class, price, y, capacity(f),'theft','booking',[]);
end

% SBLP
for f=1:num_flight
    flight_sblp_standard(f) = Flight(num_class, price, y_sblp(1:end-1,f),capacity(f),'standard','booking',[]);
    flight_sblp_theft(f) = Flight(num_class, price, y_sblp(1:end-1,f),capacity(f),'theft','booking',[]);
    flight_sblp_assortment(f) = Flight(num_class, price, zeros(num_class-1,1),capacity(f),'standard','assortment',squeeze(assortment_control(:,f,:,:)));
end

% Updated
for f = 1:num_flight
    flight_updated_theft(f) = Flight(num_class, price, zeros(num_class-1,1),capacity(f),'theft','booking',[]);
    flight_updated_standard(f) = Flight(num_class, price, zeros(num_class-1,1),capacity(f),'standard','booking',[]);
end

%
demand = zeros(num_class,num_flight,num_training); % this is for EMSR-b

%% Update standard nesting controls
Kapacity = repmat(capacity',num_class - 1,1);
Y = y_sblp(1:end-1,:);
B = Kapacity - Y;
revenue_update_standard = zeros(num_training,1);
for h=1:num_training
    h
    for f=1:num_flight
        put_protection_level(flight_updated_standard(f),Y(:,f));
        reset(flight_updated_standard(f));
    end
    [r,g, demand(:,:,h)] = evaluate_sample_path_gradient_standard_one(num_class,num_flight,num_channel,num_period,arrival_rate,v,w,v0,price,flight_updated_standard);
    %[r,g, demand(:,:,h)] = evaluate_sample_path_gradient_standard_one_fixed(num_class,num_flight,num_channel,num_period,arrival_rate,v,w,v0,price,flight_updated_standard);
    %[r,g, demand(:,:,h)] = evaluate_sample_path_gradient_standard(num_class,num_flight,num_channel,num_period,arrival_rate,V,W,V0,price,flight_updated_standard);
    % [r,g, demand(:,:,h)] = evaluate_sample_path_gradient_standard_group(num_class,num_flight,num_channel,num_period,arrival_rate,V,W,V0,price,flight_updated_standard);
    %[r,g, demand(:,:,h)] = evaluate_sample_path_gradient_standard_group_attractiveness(num_class,num_flight,num_channel,num_period,arrival_rate,v,w,v0,price,flight_updated_standard);
    revenue_update_standard(h) = r;
    B = min(B+squeeze(g)/(100*h+2000),Kapacity);
    Y = Kapacity - B;
    Y = min(Y, Kapacity);
end
figure;plot(revenue_update_standard);
xlabel('Sample Path');ylabel('Revenue');title('Updating Standard Nesting Controls')
% load Y;
for f = 1:num_flight
    put_protection_level(flight_updated_standard(f),Y(:,f));
end

%% Update theft nesting controls
Y = y_sblp(1:end-1,:);
revenue_update_theft = zeros(num_training,1);
for h=1:num_training
    h
    for f=1:num_flight
        put_protection_level(flight_updated_theft(f),Y(:,f));
        reset(flight_updated_theft(f));
    end
    [r,g, demand(:,:,h)] = evaluate_sample_path_gradient_theft_one(num_class,num_flight,num_channel,num_period,arrival_rate,v,w,v0,price,flight_updated_theft);
    %[r,g, demand(:,:,h)] = evaluate_sample_path_gradient_theft_one_fixed(num_class,num_flight,num_channel,num_period,arrival_rate,v,w,v0,price,flight_updated_theft);
    %[r,g, demand(:,:,h)] = evaluate_sample_path_gradient_theft(num_class,num_flight,num_channel,num_period,arrival_rate,V,W,V0,price,flight_updated_theft);
    % [r,g, demand(:,:,h)] = evaluate_sample_path_gradient_theft_group(num_class,num_flight,num_channel,num_period,arrival_rate,V,W,V0,price,flight_updated_theft);
    %[r,g, demand(:,:,h)] = evaluate_sample_path_gradient_theft_group_attractiveness(num_class,num_flight,num_channel,num_period,arrival_rate,v,w,v0,price,flight_updated_theft);
    revenue_update_theft(h) = r;
    Y = max(Y+squeeze(g)/(100*h+2000),0);
    Y = min(Y, Kapacity);
end
figure;plot(revenue_update_theft);
xlabel('Sample Path');ylabel('Revenue');title('Updating Theft Nesting Controls')
% load Y;
for f = 1:num_flight
    put_protection_level(flight_updated_theft(f),Y(:,f));
end

% EMSRb
for f=1:num_flight
    mu = squeeze(mean(demand(:,f,:),3));
    sigma = squeeze(std(demand(:,f,:),0,3));
    y = get_EMSRb(num_class,price, mu, sigma, capacity(f));
    flight_EMSRb_standard(f) = Flight(num_class, price, y, capacity(f),'standard','booking',[]);
    flight_EMSRb_theft(f) = Flight(num_class, price, y, capacity(f),'theft','booking',[]);
end