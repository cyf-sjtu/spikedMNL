function [r,g, demand] = evaluate_sample_path_gradient_theft_one(numClass,numFlight,numChannel,numEpoch,rate,v,w,v0,fare,flight)
    demand = zeros(numClass,numFlight);
    num = 1e4;
    R = zeros(num+1,1); % revenue
    numProduct = numClass*numFlight;
    % seq = zeros(num+1,2); % product chosen [j f]
    pRpx = zeros(numFlight,num+1); % partial R partial x, marginal value of inventory
    pRpy = zeros(numProduct,num+1); % parial R partial y, marginal value of control
    % x = sparse(zeros(numFlight,num+1)); % inventory
    u = zeros(numClass,numFlight,num+1); % accepted request
    pupx = zeros(numFlight,numProduct,num+1); % partial u partial x
    pupy = zeros(numProduct,numProduct,num+1); % partial u partial y
    % forward
    tau = 0;
    for t=1:numEpoch
        seqArrival = generate_multichannel_arrival(numChannel,rate(:,t),1.0);
        v_ct = v(:,:,:,t);
        v_ct = reshape(v_ct,numProduct,numChannel);
        w_ct = w(:,:,:,t);
        w_ct = reshape(w_ct,numProduct,numChannel);
        for n=1:size(seqArrival,1)
            c = seqArrival(n,2);
            % seed = randi(10^9); 
            % pseudo flights
            q = 1.0;
            I = [];
            [avail,spike] = get_avail_spike(numFlight,numClass,flight,c,t);
            vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
            p=choose_product(numFlight*numClass,v0(c,t),vtemp,avail);
            if(p > 0)
                tau = tau + 1;
                demand(p) = demand(p) + q;
%             end
            while(q > 1e-9)
                if(p > 0)
                    I = [I,p];
                    [f,j] = find_flight_product(p,numFlight,numClass);
                    [u(j,f,tau), ind1, ind2] = put_pseudo_booking(flight(f),j,q);
                    q = q - u(j,f,tau);
                    if (ind2 == 1)
                        for ii = I(1:end-1)
                            %[ff,jj] = find_flight_product(ii,numFlight,numClass);
                            pupx(:,p,tau) = pupx(:,p,tau) - pupx(:,ii,tau);
                            pupy(:,p,tau) = pupy(:,p,tau) - pupy(:,ii,tau);
                        end
                        break;
                    elseif (ind1 == 2)
                        pupx(f,p,tau) = 1;
                        if(j~=1)
                            pupy(p,p,tau) = -1;
                        end
                        for ii= I(1:end-1)
                            [ff,jj] = find_flight_product(ii,numFlight,numClass);
                            if(f==ff)
                                pupx(:,p,tau) = pupx(:,p,tau) - pupx(:,ii,tau);
                                pupy(:,p,tau) = pupy(:,p,tau) - pupy(:,ii,tau);
                            end
                        end
                    end
                else
                    break;
                end
                [avail,spike] = get_avail_spike(numFlight,numClass,flight,c,t);
                vtemp = v_ct(:,c); vtemp(logical(spike)) = w_ct(logical(spike),c);
                p=choose_product(numFlight*numClass,v0(c,t),vtemp,avail);
            end
            end
        end
    end
    % backward
    price1 = repmat(fare',numFlight,numFlight);
    price2 = repmat(fare',numProduct,numFlight);
    for t=tau:-1:1
        R(t) = sum(u(:,:,t),2)'*fare+R(t+1);
        temp=repmat(pRpx(:,t+1)',numClass,1);
        temp=temp(:);
        pRpx(:,t) = sum(price1.*pupx(:,:,t),2)...
                    - sum(pupx(:,:,t).*repmat(temp',numFlight,1),2)...
                    + pRpx(:,t+1);
        pRpy(:,t) = sum(price2.*pupy(:,:,t),2)...
                    -sum(pupy(:,:,t).*repmat(temp',numProduct,1),2) ...
                    + pRpy(:,t+1);  
    end
    temp = reshape(pRpy(:,1),numClass,numFlight);
    g = temp(2:end,:);
    r = R(1); 
end