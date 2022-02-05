function [avail,spike] = get_avail_spike(numFlight,numClass,flight,channel,epoch)
avail = zeros(numFlight*numClass,1);
spike = zeros(numFlight*numClass,1);
for f=1:numFlight
    temp_avail = get_avail(flight(f),channel,epoch);
    avail((f-1)*numClass+1:f*numClass) = temp_avail;
    temp=find(temp_avail == 1); 
    if(~isempty(temp) )
        idx = temp(end); 
        spike((f-1)*numClass+idx) = 1;
    end
end
end