function y = get_robust_protection_level(m,fare,K)
    y = [];
    if(m==length(fare))
        delta = m - sum(fare(2:end)./fare(1:end-1));
        y = zeros(m-1,1);
        for i=1:m-1
            y(i) = K/delta*(i - sum(fare(2:i+1)./fare(1:i)));
        end
    end
end
