classdef planet
    properties
        radius;
        mass;
        eccentricity;
        volume;
        star;
        nplanets;
        dist2star;
    end
    methods
        function obj = planet(r,m,s)
            obj.radius = r;
            obj.mass = m;
            obj.star = s;
        end
    end
end
