function seqArrival = generate_multichannel_arrival(numChannel,rate,lengEpoch)
    numArrival=poissrnd(rate); % not the actual arrival
    seqArrival = zeros(sum(numArrival),2);
    pos=1;
    for c=1:numChannel
        seqArrival(pos:pos+numArrival(c)-1,1) = rand(numArrival(c),1)*lengEpoch; % arrival time
        seqArrival(pos:pos+numArrival(c)-1,2) = c;
        pos = pos+numArrival(c);
    end
    seqArrival = sortrows(seqArrival);
end
