function [U, I] = generate_customer(V,V0,seed)
    rng(seed);
    e = -log(-log(rand(length(V)+1,1)))-0.5772;
    U = [V;V0]+e;
    [B,I] = sort(U,'descend');
    len = find(I==length(V)+1);
    I = I(1:len-1);
    U = U(1:len-1);
end