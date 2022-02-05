function [k,j] = find_flight_product(index,f,m)
    k = floor((index-0.1)/m) + 1; % flight index
    if(k > f) 
        k=0;
        j=0; % no purchase
    end
    j = mod(index,m) ; % product index
    if (j==0)
        j=m;
    end
end