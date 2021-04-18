classdef commands
    methods (Static)
        
        %Longitud del arco de meridiano----------------------------
        function lamlat = lat2lam(a,invf,lat)
            %lat2lam Calcula la latitud de arco de meridiano a partir de la
            %latitud geodesica
            lat = lat*pi/180;
            [f,e2,e] = commands.cteElli(a,invf);
            [g1,g2,g3,g4,g5,g6] = commands.g(e2,lat);
            lamlat = a*(1-e2)*(g1+g2+g3+g4+g5+g6);
            
        end
        
        function lat = lam2lat(a,invf,lam)
            %lam2lat calcula la latitud geodesica a partir de la latitud de
            %arco de meridiano
            
            [f,e2,e] = commands.cteElli(a,invf);
            lat0 = lam/(a*(1-e2));
            latdif = 1;
            threshold = 0.000001/60/60;
            while latdif>threshold
                [g1,g2,g3,g4,g5,g6] = commands.g(e2,lat0);
                lat1 = lam/(a*(1-e2))-(g2+g3+g4+g5+g6);
                latdif = abs(lat1-lat0)/pi*180;
                lat0 = lat1;
                
            end
            lat = lat1*180/pi;
            
        end
        
        
        %Latitud creciente----------------------------------------
        
        function q = lat2q(a,invf,lat)
            %lat2q Calcula la latitud creciente a partir de la latitud
            %geodesica
            lat = lat*pi/180;
            [f,e2,e] = commands.cteElli(a,invf);
            q = log(((1-e*sin(lat))/(1+e*sin(lat)))^(e/2)*tan(lat/2+pi/4));
        end
        
        function lat = q2lat(a,invf,q)
            %q2lat Calcula la latitud geodesica a partir de la latitud
            %creciente
           [f,e2,e] = commands.cteElli(a,invf);
           lat0 = 2*(atan(exp(q))-pi/4);
           threshold = (0.000001/60/60)*pi/180;
           latdif = 1; 
           while latdif > threshold
                lat1 = 2*(atan(exp(q)*((1+e*sin(lat0))/(1-e*sin(lat0)))^(e/2))-pi/4);
                latdif = abs(lat1-lat0);
                lat0 = lat1;
           end
           lat = lat1*180/pi;
        end
        
        
        %Latitud autalica-----------------------------------------
        
        function aut = lat2aut(a,invf,lat)
            %lat2aut Calcula la latitud autalica a partir de la latitud
            %geodesica
            lat = lat*pi/180;
            [f,e2,e] = commands.cteElli(a,invf);
            aut = ((1-e2)/(4*e))*((2*e*sin(lat))/(1-e2*(sin(lat))^2)+log((1+e*sin(lat))/(1-e*sin(lat))));
        end
        
        function lat = aut2lat(a,invf,aut)
            %aut2lat Calcula la latitud geodesica a partir de la latitud
            %autalica
           [f,e2,e] = commands.cteElli(a,invf);
           lat0 = asin(aut);
           threshold = (0.000001/60/60)*pi/180;
           latdif = 1; 
           while latdif > threshold
               f = (1-e2)/(4*e)*((2*e*sin(lat0))/(1-e2*sin(lat0)^2)+log((1+e*sin(lat0))/(1-e*sin(lat0))))-aut;
               f2 = ((1-e2)*cos(lat0))/((1-e2*sin(lat0)^2)^2);
               lat1 = lat0-f/f2;
               latdif = abs(lat1-lat0);
               lat0 = lat1; 
           end
           lat = lat1*180/pi;
        end 
    end
    
    
    %metodos privados auxiliares..........................................
    
    methods (Access = private,Static)
        function [g1,g2,g3,g4,g5,g6] = g(e2,lat)
            e = sqrt(e2);
            g1 = lat;
            g2 = (3*e^2*(lat-sin(lat)*cos(lat)))/4;
            g3 = (15/256)*e^4*(12*lat-8*sin(2*lat)+sin(4*lat));
            g4 = (35/3072)*e^6*(60*lat-45*sin(2*lat)+9*sin(4*lat)-sin(6*lat));
            g5 = (315/393216)*e^8*(840*lat-672*sin(2*lat)+168*sin(4*lat)-32*sin(6*lat)+3*sin(8*lat));
            g6 = (693/2621440)*e^10*(2520*lat-2100*sin(2*lat)+600*sin(4*lat)-150*sin(6*lat)+25*sin(8*lat)-2*sin(10*lat));
        end
        
        function [g1,g2,g3,g4,g5,g6] = gaut(e2,lat)
            e = sqrt(e2);
            g1 = 1;
            g2 = 2/3*e^2*sin(lat)^2;
            g3 = 3/5*e^4*sin(lat)^4;
            g4 = 4/7*e^6*sin(lat)^6;
            g5 = 5/9*e^8*sin(lat)^8;
            g6 = 6/11*e^10*sin(lat)^10;
        end
        
        function [f,e2,e] = cteElli(a,invf)
            f = 1/invf;
            b = a*(1-f);
            e2 = (a^2-b^2)/a^2;
            e = sqrt(e2);
        end
    end
end