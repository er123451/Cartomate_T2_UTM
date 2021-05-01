classdef commands2
    %COMMANDS2 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static)
        function [x,y] = GKdir(long,lat,long0,k0,tx,ty,a,invf)
            %GKdir problema directo de la proyección de Gauß-Krüger
            %mediante el método aproximado.
            long = long*pi/180;
            lat = lat*pi/180;
            long0 = long0*pi/180;
            Along = long-long0;
            [q0,q1,q2,q3,q4,q5,q6,q7,q8,q9] = commands2.fq(lat,a,invf);       
            x = tx +k0*(q1*Along-q3*Along^3/6+q5*Along^5/120-q7*Along^7/5040+q9*Along^9/362880);
            y = ty +k0*(q0-q2*Along^2/2+q4*Along^4/24-q6*Along^6/720+q8*Along^8/40320);
        end
        
        function [lat,long] = GKinv(x,y,long0,k0,tx,ty,a,invf)
           x = (x-tx)/k0;
           y = (y-ty)/k0;
           
        end
        
    end
    methods (Access = private, Static)
        function [f1,h2,f3,h4,f5,h6,f7,h8,f9] = fy(lat,a,invf)
           %-----------cte--------------------------
           [f,e2,e] = commands.cteElli(a,invf);
           ei = e/sqrt(1-e2);
           lat = commands.lam2lat(a,invf,y);
           v =  a*(1-e2*(sin(lat))^2)^(-1/2);
           t = tan(lat);
           n = ei*cos(lat);
           %---------------valores y----------------
           f1 = (v*cos(lat))^-1;
           h2 = t*v^-2*(1+n^2);
           f3 = (v^3*cos(lat))^-1*(1+2*t^2+n^2);
           h4 = t*v^-4*(5+3*t^2+n^2*(6-3*n^2)-6*n^2+t^2);
           f5 = (v^5*cos(lat))^-1*(5+t^2*(28+24*t^2)+6*n^2+8*n^2*t^2);
           h6 = t*v^-6*(61+t^2*(90+45*t^2));
           f
           
        end
        function [q0,q1,q2,q3,q4,q5,q6,q7,q8,q9] = fq(lat,a,invf)
           %-----------cte--------------------------
           [f,e2,e] = commands.cteElli(a,invf);
           ei = e/sqrt(1-e2);
           v =  a*(1-e2*(sin(lat))^2)^(-1/2);
           t = tan(lat);
           n = ei*cos(lat);
           %---------------valores q----------------
           q0 = commands.lat2lam(a,invf,lat*180/pi);
           q1 = v*cos(lat);
           q2 = -v*t*cos(lat)^2;
           q3 = -v*cos(lat)^3*(1-t^2+n^2);
           q4 = v*t*cos(lat)^4*(5-t^2+n^2*(9+4*n^2));
           q5 = v*cos(lat)^5*(5-t^2*(18-t^2)+n^2*(14+13*n^2)-n^2*t^2*(58+64*n^2));
           q6 = -v*t*cos(lat)^6*(61-t^2*(58-t^2)+n^2*(270+445*n^2)-330*n^2*t^2);
           q7 = -v*cos(lat)^7*(61-t^2*(479-179*t^2+t^4)+331*n^2-3298*n^2*t^2+1771*n^2*t^4);
           q8 = v*t*cos(lat)^8*(1385-t^2*(3111-543*t^2)+10899*n^2-32802*n^2*t^2);
           q9 = v*cos(lat)^9*(1385-t^2*(19028-18270*t^2));
        end
    end
end

