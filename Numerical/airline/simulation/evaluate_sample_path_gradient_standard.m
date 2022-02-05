function [r,g, demand] = evaluate_sample_path_gradient_standard(numClass,numFlight,numChannel,numEpoch,rate,V,W,V0,fare,flight)
    demand = zeros(numClass,numFlight);
    num = 1e4;
    R = zeros(num+1,1); % revenue
    % seq = zeros(num+1,2); % product chosen [j f]
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
        V_ct = V(:,:,:,t);
        V_ct = reshape(V_ct,numClass*numFlight,numChannel);
        W_ct = W(:,:,:,t);
        W_ct = reshape(W_ct,numClass*numFlight,numChannel);
        for n=1:size(seqArrival,1)
            c = seqArrival(n,2);
            seed = randi(10^9); 
            % pseudo flights
            [avail,spike] = get_avail_spike(numFlight,numClass,flight,c,t);
            Vtemp = V_ct(:,c); Vtemp(logical(spike)) = W_ct(logical(spike),c);
            % [U,q,I] = generate_group_customer(Vtemp,V0(c,t),avail,rate(c,t),seed);
            [U, I] = generate_customer(Vtemp,V0(c,t),seed);
            q = 1.0;
            if(~isempty(I))
                demand(I(1)) = demand(I(1)) + q;
                tau=tau+1;
                for i = 1:length(I)
                    [f,j] = find_flight_product(I(i),numFlight,numClass);
                    [u(j,f,tau), ind1, ind2] = put_pseudo_booking(flight(f),j,q);
                    q = q - u(j,f,tau);
                    if (ind2 == 1)
                        for ii=1:i-1
                            [ff,jj] = find_flight_product(I(ii),numFlight,numClass);
                            %pupx(:,j,f,tau) = pupx(:,j,f,tau) - pupx(:,jj,ff,tau);
                            pupb(:,:,j,f,tau) = pupb(:,:,j,f,tau) - pupb(:,:,jj,ff,tau);
                            pups(:,:,j,f,tau) = pups(:,:,j,f,tau) - pups(:,:,jj,ff,tau);
                        end
                        break;
%                     elseif (ind2 == 3)
%                         pupx(f,j,f,tau) = 1;
%                         for ii=1:i-1
%                             [ff,jj] = find_flight_product(I(ii),numFlight,numClass);
%                             if(f==ff)
%                                 pupx(:,j,f,tau) = pupx(:,j,f,tau) - pupx(:,jj,ff,tau);
%                                 pupb(:,:,j,f,tau) = pupb(:,:,j,f,tau) - pupb(:,:,jj,ff,tau);
%                                 pups(:,:,j,f,tau) = pups(:,:,j,f,tau) - pups(:,:,jj,ff,tau);
%                             end
%                         end
                    elseif (ind1 == 2)
                        %if(j~=1)
                            pupb(j,f,j,f,tau) = 1;
                        %end
                        pups(j:end,f,j,f,tau) = -1;
                        for ii=1:i-1
                            [ff,jj] = find_flight_product(I(ii),numFlight,numClass);
                            if(f==ff && jj>j)
                                %pupx(:,j,f,tau) = pupx(:,j,f,tau) - pupx(:,jj,ff,tau);
                                pupb(:,:,j,f,tau) = pupb(:,:,j,f,tau) - pupb(:,:,jj,ff,tau);
                                pups(:,:,j,f,tau) = pups(:,:,j,f,tau) - pups(:,:,jj,ff,tau);
                            end
                        end
                    end
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
    g = pRpb(2:end,:,1);
    r = R(1); 
end