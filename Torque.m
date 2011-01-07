% Evaluating The Torque between 2 parallel PM's
% allan liu 7/01/2011

function [Tx,Ty,Tz]= Torque(a1,b1,c1,a2,b2,c2,a,b,c,d,e,f,b1,b2);  

%a1...c1 are the half lengths of pm1 in meters
%a2...c2 are the half lengths of pm2
%a,b,c is the coordinates of pm2 relative to pm1 in meters
%d,e,f is the coordinates of pm1 relative to Ot
%b1 and b2 are the flux density components along the z axis in tesla

Tx=0;   %initialising Torque vector
Ty=0;
Tz=0; 
mo=4*pi*10^-7 %constant

for ii=[0,1]
    for jj=[0,1]
        for kk=[0,1]
            for ll=[0,1]
                for mm=[0,1]
                    for nn=[0,1]
                        
                        
                        
                        Ex=
                        Tx=Tx+((b1*b2/(4*pi*mo))*(-1).^(ii+jj+kk+ll+mm+nn)*Ex)
                        
                      



end

