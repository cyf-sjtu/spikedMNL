function y = get_EMSRb(m,fare, mu, sigma, K)
    y=[];
    if(m==length(fare) && m==length(mu) && m==length(sigma))
        y=zeros(m-1,1);
        weighted_fare = zeros(m-1,1);
        for j=1:m-1
            weighted_fare(j) = sum(fare(1:j).*mu(1:j))/sum(mu(1:j));
            y(j) = min([norminv(1-fare(j+1)/weighted_fare(j), sum(mu(1:j)), sqrt(sum(sigma(1:j).^2)));K]);
        end
    end
end