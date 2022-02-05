classdef Flight < handle
    properties(SetAccess = private)
        % array
        numClass;   % # of fare classes
        protectionLevel;  % protection level
        bookingLimit; % booking limit
        strategy; % 'standard' or 'theft' nesting
        policy; % 'booking' or 'assortment'
        timeControl; %
        capacity;
        inventory; 
        bookings;
        cumulativeBookings;
        reverseCumulativeBookings;
        fare;
        revenue;
        iscell;
    end
    methods
        function obj = Flight(m, fare, y, K, s, p, t)
            if nargin > 0
                % initialize
                obj.numClass = m;
                obj.fare = fare;
                obj.protectionLevel = [0;y]; % shifted
                obj.bookingLimit = K-[0;y];
                obj.strategy = s;
                obj.capacity = K;
                obj.inventory = obj.capacity;
                obj.bookings = zeros(m,1);
                obj.cumulativeBookings = zeros(m,1);
                obj.reverseCumulativeBookings = zeros(m,1);
                obj.revenue = 0;
                obj.policy = p;
                obj.timeControl = t;
                obj.iscell = iscell(t);
            end
        end
        function reset(obj)
            obj.inventory = obj.capacity;
            obj.bookings = zeros(obj.numClass,1);
            obj.cumulativeBookings = zeros(obj.numClass,1);
            obj.reverseCumulativeBookings = zeros(obj.numClass,1);
            obj.revenue = 0;
        end
        function put_protection_level(obj,y)
            obj.protectionLevel = [0;y]; % shifted
            obj.bookingLimit = obj.capacity - obj.protectionLevel;
        end
        function avail = get_avail(obj,c,t)
            avail = zeros(obj.numClass,1);
            if (obj.inventory<=0) 
                return;
            end
            switch obj.policy
                case 'booking'
                switch obj.strategy
                    case 'standard'
                        avail= (obj.reverseCumulativeBookings < obj.bookingLimit);
                    case 'theft'
                        avail = (obj.inventory > obj.protectionLevel);
                    otherwise
                        warning('Unexpected strategy!')
                end
                case 'assortment'
                    if(obj.iscell)
                        time = rand();
                        for j = 1:obj.numClass
                            len = size(obj.timeControl{j,c,t},1);
                            for l = 1:len
                                avail(j) = avail(j)|(time>= obj.timeControl{j,c,t}(l,1) && time<= obj.timeControl{j,c,t}(l,2));
                            end
                        end 
                    else
                        avail = (obj.timeControl(:,c,t) >= rand(obj.numClass,1));
                    end
            end
        end
        function put_booking(obj,j)
            switch obj.policy
                case 'booking'
                switch obj.strategy
                    case 'standard'
                        check = (obj.inventory>0) && (obj.reverseCumulativeBookings(j) < obj.bookingLimit(j));
                    case 'theft'
                        check = (obj.inventory > obj.protectionLevel(j));
                    otherwise
                        warning('Unexpected strategy!')
                end
                case 'assortment'
                    check = 1;
            end
            if(check)
                obj.inventory = obj.inventory - 1;
                obj.bookings(j) = obj.bookings(j)+1;
                obj.cumulativeBookings = cumsum(obj.bookings);
                obj.reverseCumulativeBookings = cumsum(obj.bookings,'reverse');
                obj.revenue = obj.revenue + obj.fare(j);
            end
        end
        % for sample approximation
        function [u, ind1, ind2] = put_pseudo_booking(obj,j,q)
            u = 0; ind1 = 0; ind2 = 0;
            switch obj.policy
                case 'booking'
                switch obj.strategy
                    case 'standard'
                        check = (obj.inventory>0) && (obj.reverseCumulativeBookings(j) < obj.bookingLimit(j));
                    case 'theft'
                        check = (obj.inventory > obj.protectionLevel(j));
                    otherwise
                        warning('Unexpected strategy!')
                end
                case 'assortment'
                    check = 1;
            end
            if(check)
                switch obj.strategy
                    case 'standard'
                        [temp, ind1] = max([0.0,obj.bookingLimit(j) - obj.reverseCumulativeBookings(j)]);
                        [u, ind2] = min([q,obj.inventory,temp]);
                    case 'theft'
                        [temp, ind1] = max([0.0,obj.inventory-obj.protectionLevel(j)]);
                        [u, ind2] = min([q,temp]);
                end
                obj.inventory = obj.inventory - u;
                obj.bookings(j) = obj.bookings(j) + u;
                obj.cumulativeBookings = cumsum(obj.bookings);
                obj.reverseCumulativeBookings = cumsum(obj.bookings,'reverse');
                obj.revenue = obj.revenue + u*obj.fare(j);
            end
        end
        
        function cum = get_cumulative_booking(obj)
            cum = obj.cumulativeBookings;
        end
        function level = get_protection_level(obj)
            level = obj.protectionLevel;
        end
        function revenue = get_revenue(obj)
            revenue = obj.revenue;
        end
        function sale = get_sale(obj)
            sale = obj.bookings;
        end
    end
end
