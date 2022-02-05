function [U,q,I] = generate_group_customer(V,V0,avail,rate,seed)
    rng(seed);
    e = -log(-log(rand(length(V)+1,1)))-0.5772;
    U = [V;V0]+e;
    [B,I] = sort(U,'descend');
    len = find(I==length(V)+1);
    I = I(1:len-1);
    U = U(1:len-1,:);
    v = exp(V)'*avail;
    v0 = exp(V0);
    q = v/(v+v0) * poissrnd(rate);
    % N = N(I);
end