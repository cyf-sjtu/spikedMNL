function p=choose_product(numProduct,V0,V,availability)
temp = 0;
total_v = V0+availability'*V;
product=1:numProduct;
product=product(logical(availability));
r=rand();
p=-1;   % default
for j=product
    temp = temp + V(j)/total_v;
    if temp >= r
        p = j; break;
    end
end
end