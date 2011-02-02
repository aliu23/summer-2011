% Evaluating The Torque between 2 parallel PM's
% allan liu 7/01/2011 

function [Tx,Ty,Tz]= Torque(a1,b1,c1,a2,b2,c2,a,b,c,d,e,f,br1,br2)  

%a1...c1 are the half lengths of pm1 in meters
%a2...c2 are the half lengths of pm2
%a,b,c is the coordinates of pm2 relative to pm1 in meters
%d,e,f is the coordinates of pm1 relative to Ot
%b1 and b2 are the flux density components along the z axis in tesla

Tx=0;   %initialising Torque vector
Ty=0;
Tz=0; 
mo=4*pi*10^-7; %constant

for ii=[0,1]
  for jj=[0,1]
    for kk=[0,1]
      for ll=[0,1]
        for mm=[0,1]
          for nn=[0,1]
            
            Cw=(-1)^mm.*c1-f;
            Cv=(-1)^kk.*b1-e;
            Cu=(-1)^ii.*a1-d;
            
            w=c-(-1)^mm.*c1+(-1)^nn.*c2;
            v=b-(-1)^kk.*b1+(-1)^ll.*b2;
            u=a-(-1)^ii.*a1+(-1)^jj.*a2;
            
            s=sqrt(u.^2+v.^2+w.^2);
            
            Ez=(1/36).*(-u.^3-18.*v.*u.^2-6.*u.*(w.^2+6.*Cu...
              .*v-3.*v.*(2.*Cv+v)+3.*Cv.*s)+v.*(v.^2+6.*(w.^2+...
              3.*Cu.*s))+6.*w.*(w.^2-3.*v.*(2.*Cv+v)).*atan(u./w)...
              -6.*w.*(w.^2-3.*u.*(2.*Cu+u)).*atan(v./w)-9.*...
              (2.*(v.^2+2.*Cv.*v-u.*(2.*Cu+u)).*w.*atan(u.*v./(w.*s))...
              -2.*u.*(2.*Cu+u).*v.*log(s-u)-(2.*Cv+v).*(v.^2-w.^2)...
              .*log(u+s)+2.*u.*v.*(2.*Cv+v).*log(s-v)+(2.*Cu+...
              u).*(u.^2-w.^2).*log(v+s)));
            
            Ey=(1/8)*...
              ((2.*Cw+w).*u.^2-8.*u.*v.*(Cw+w)+8.*u.*v.*(Cw+w).*log(s-v)...
              +4.*Cw.*u.*s+6.*w.*s.*u+(2.*Cw+w).*(v.^2+w.^2)+...
              4.*w.*(w.^2+2.*Cw.*w-u.*(2.*Cu+u)).*acoth(u./s)+...
              4.*v.*(-2.*Cu.*w.*acoth(v./s)+2.*w.*(Cw+w).*atan(u./w)...
              +(w.^2+2.*Cw.*w-u.*(2.*Cu+u)).*atan(u.*v./(w.*s)))...
              -2.*(2.*Cw+w).*(v.^2+w.^2).*log(u+s)+8.*Cu.*w.*s);
            
            Ex=(1/8).*(-2.*Cw.*(-4.*v.*u+s.^2+2.*v.*s)-w.*...
              (-8.*v.*u+s.^2+8.*Cv.*s+6.*v.*s)+4.*(2.*Cv.*u.*...
              w.*acoth(u./s)+w.*(v.^2+2.*Cv.*v-w.*(2.*Cw+w))...
              *acoth(v./s)-u.*(2.*w.*(Cw+w).*atan(v./w)+2*v.*...
              (Cw+w).*log(s-u)+(w.^2+2.*Cw.*w-v.*(2.*Cv+v)).*...
              atan(u.*v./(w.*s))))+2.*(2.*Cw+w).*(u.^2+w.^2).*log(v+s));
            
            Tx=Tx+(-1)^(ii+jj+kk+ll+mm+nn)*Ex;
            Ty=Ty+(-1)^(ii+jj+kk+ll+mm+nn)*Ey;
            Tz=Tz+(-1)^(ii+jj+kk+ll+mm+nn)*Ez;
            
          end
        end
      end
    end
  end
end

Tx=Tx.*br1*br2/(4*pi*mo);
Ty=Ty.*br1*br2/(4*pi*mo);
Tz=Tz.*br1*br2/(4*pi*mo);
                        
end
