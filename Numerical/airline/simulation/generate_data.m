% rng(20171019)
num_period = 100;
num_channel = 5;
arrival_rate = 0.5 + 1.5 * rand(num_channel, num_period); % num_channel-by-num_period  
num_class = 4;
num_flight = 10;
capacity = 25*ones(num_flight,1);
price = sort(randi([50,500],num_class,1),'descend');

v = 10*rand(num_class,num_flight,num_channel,num_period);
v = sort(v);
w = v + 2.0*rand(num_class,num_flight,num_channel,num_period);
% w = v;

v0 = squeeze(sum(sum(v))).* (0.7+0.3*sort(rand(num_channel,num_period),2,'descend'));

V = log(v); W = log(w); V0=log(v0);