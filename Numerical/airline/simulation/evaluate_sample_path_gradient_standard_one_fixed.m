function [r,g, demand] = evaluate_sample_path_gradient_standard_one_fixed(numClass,numFlight,numChannel,numEpoch,rate,v,w,v0,fare,flight)
    demand = zeros(numClass,numFlight);
    num = 1e4;
    R = zeros(num+1,1); % revenue
    %seq = zeros(num+1,2); % product chosen [j f]
    % pRpx = zeros(numFlight,num+1); % partial R partial x, marginal value of inventory
    pRpb = zeros(numClass,numFlight,num+1); % parial R partial b, marginal value of control
    pRps = zeros(numClass,numFlight,num+1); % parial R partial s, marginal value of sales
    % x = sparse(zeros(numFlight,num+1)); % inventory
    u = zeros(numClass,numFlight,num+1); % accepted request
    % pupx = zeros(numFlight,numClass,numFlight,num+1); % partial u partial x
    pupb = zeros(numClass,numFlight,numClass,numFlight,num+1); % partial u partial b
    pups = zeros(numClass,numFlight,numClass,numFlight,num+1); % partial u partial s
    % forward
    tau = 0;
    for t=1:numEpoch
        seqArrival = generate_multichannel_arrival(numChannel,rate(:,t),1.0);
        v_ct = v(:,:,:,t);
        v_ct = reshape(v_ct,numClass*numFlight,numChannel);
        w_ct = w(:,:,:,t);
        w_ct = reshape(w_ct,numClass*numFlight,numChannel);
        for n=1:size(seqArrival,1)
            c = seqArrival(n,2);
            % seed = randi(10^9); 
            % pseudo flights
            q = 1.0;
            %I = [];
            [avail,spike] = get_avail_spike(numFlight,numClass,flight,c,t);
            vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
            p = choose_product(numFlight*numClass,v0(c,t),vtemp,avail);
            if(p > 0)
                tau = tau + 1;
                demand(p) = demand(p) + q;
                [f,j] = find_flight_product(p,numFlight,numClass);
                %seq(tau,1) = j;seq(tau,2) = f;
                [u(j,f,tau), ind1, ind2] = put_pseudo_booking(flight(f),j,q);
                %q = q - u(j,f,tau);
                if (ind2 == 1)
%                 elseif (ind2 == 2)
%                     pupx(f,j,f,tau) = 1;
                elseif (ind1 == 2)
                    pupb(j,f,j,f,tau) = 1;
                    pups(j:end,f,j,f,tau) = -1;
                end
            end
        end
    end
    % backward
    for t=tau:-1:1
        R(t) = sum(u(:,:,t),2)'*fare+R(t+1);
%         pRpx(:,t) = pRpx(:,t+1);
%         pRpb(:,:,t) = pRpb(:,:,t+1);
%         pRps(:,:,t) = pRps(:,:,t+1);
        %j = seq(t,1);
        %f = seq(t,2);
        for f=1:numFlight
%             pRpx(f,t) = sum(sum(repmat(fare,1,numFlight).*squeeze(pupx(f,:,:,t))))...
%                     - sum(squeeze(pupx(f,:,:,t)))*pRpx(:,t+1)+ pRpx(f,t+1) ...
%                     + sum(sum(squeeze(pupx(f,:,:,t)).*pRps(:,:,t+1)));
            for j=1:numClass
                pRpb(j,f,t) =sum(sum(repmat(fare,1,numFlight).*squeeze(pupb(j,f,:,:,t))))...
                    + pRpb(j,f,t+1)...
                    + sum(sum(squeeze(pupb(j,f,:,:,t)).*pRps(:,:,t+1)));
                pRps(j,f,t) =sum(sum(repmat(fare,1,numFlight).*squeeze(pups(j,f,:,:,t))))...
                    + sum(sum(squeeze(pups(j,f,:,:,t)).*pRps(:,:,t+1)))+pRps(j,f,t+1);
            end
        end
    end
    %pRpb(1,:,1)
    g = pRpb(2:end,:,1);
    r = R(1); 
end