
function [F,E,P]=elliptic123(a1,a2,a3)
%ELLIPTIC123 computes the first, second and third elliptic integrals for
% both the complete and incomplete cases and no restriction on the input
% arguments. (Modulo some bugs; see below.)
%
% [F,E]=elliptic123(b,m)
%   Calculates incomplete elliptic integrals of the first and second kind,
%   F and E, respectively.
%    - Phase angle b may be any real or complex number
%    - Parameter m can be any real number
%
% [K,E]=elliptic123(m)
%   Calculate complete elliptic integrals of the first and second kind,
%   K and E, respectively.
%    - Equivalent to  [K,E]=elliptic12(pi/2,m)  but calculated more
%      efficiently 
%
% [K,E,P]=elliptic123(m,n)
%  Calculates the complete elliptic integrals of the first, second and
%  third kind, K, E and P respectively
%  The n parameter is only used for the elliptic PI case (third case)
%
% [F,E,P]=elliptic123(b,m,n)
%  Calculates the incomplete elliptic integrals of the first, second and
%  third kind, F, E and P respectively
%  The n parameter is only used for the elliptic PI case (third case)
%
% There is a bug in the incomplete case for F and E: 
%   when complex numbers are expected in the output, the complex part will
%   not be calculated correctly when m>(1/sin(b))^2.
% 
% There is a bug in the incomplete case for Pi 
%   when complex numbers are expected in the input, Cannot evaluate
% For the elliptic PI case when phase is greater than PI/2 and n>1, Cannot
% evaluate

if nargout<3
  
  if nargin==1
    [F,E] = elliptic12x(a1);  % == elliptic12(m) 
  elseif nargin==2
    [F,E] = elliptic12x(a1,a2); % == elliptic12(b,m)
  else
    error('Wrong number of input arguments')
  end
  
elseif nargout==3
  
  if nargin==2
    [F,E] = elliptic12x(a1); % == elliptic12(m)
    P=elliptic3x(a1,a2);   % == elliptic3(m,n)
  elseif nargin==3
    [F,E] = elliptic12x(a1,a2); % == elliptic12(b,m)
    P=elliptic3x(a1,a2,a3);   % == elliptic3(b,m,n)
  else
    error('Wrong number of input arguments')
  end
  
else
  error('Wrong number of output arguments')
end 

end

function [F,E]=elliptic12x(b,m)
%ELLIPTIC12x computes the first and second Elliptic integrals for both the
% complete and incomplete cases and no restriction on the input arguments. 
%
% [F,E]=elliptic12(b,m)
%   Calculate incomplete elliptic integrals of the first and second kind,
%   F and E, respectively.
%    - Phase angle b may be any real or complex number
%    - Parameter m can be any real number
%
% [K,E]=elliptic12(m)
%   Calculate complete elliptic integrals of the first and second kind,
%   K and E, respectively.
%    - Equivalent to  [K,E]=elliptic12(pi/2,m)  but calculated more
%      efficiently
%
% There is a bug in the incomplete case: 
%   when complex numbers are expected in the output, the complex part will
%   not be calculated correctly when m>(1/sin(b))^2.

if nargin==2
    if length(b)==1, b=b(ones(size(m))); end
    if length(m)==1, m=m(ones(size(b))); end
    if ~isequal(size(b),size(m))
        error('m and b must be equal sizes')
    end
end

%F function

 
if nargin==1  %now b= pi/2
   
    m=b;
    
    N=size(m);
    
    F=nan(size(m));
    
    if any(m<0) 
    
        mm=m(m<0);
        F(m<0)=(1./sqrt(1-mm)).*ellipke(-mm./(1-mm)); %in future may want to consider complex inputs
    
    end

    if any(m>1)
        
        mm=m(m>1);
        F(m>1)=(1./sqrt(mm)).*(elliptic12i(asin(sqrt(mm)),1./mm)); 
        
    end

    if any(m<=1&m>=0) %0<m<=1
        
        mm=m(m<=1&m>=0);
        F(m<=1&m>=0)=ellipke(mm); 
        
    end
    

elseif nargin==2 %2 inputs of b and m
    

    phase_ind = b>pi/2 | b<0; 
    mone_ind= m==1;
    
    if any(phase_ind & ~mone_ind)
        
        mm = m(phase_ind);
        bb = b(phase_ind);
        
        phi = mod(bb+pi/2,pi)-pi/2;
        a = round(bb./pi);
        F(phase_ind) = 2.*a.*elliptic12x(mm) + sign(phi).*elliptic12x(abs(phi),mm);
    
    end

    if any(mone_ind)
        F(mone_ind)=inf;
    end 
    
    mneg_ind = m<0 & ~phase_ind; %m<0 and when b IS in normal range
    
    if any(mneg_ind)
        
        mm=m(mneg_ind);
        bb=b(mneg_ind);
    
        t=asin((sin(bb).*sqrt(1-mm))./sqrt(1-mm.*(sin(bb)).^2)); %t is actually theta
        F(mneg_ind)=(1./sqrt(1-mm)).*elliptic12i(t,-mm./(1-mm)); 
    
    end
    
    mpos_ind = m>1 & ~phase_ind;
    if any(mpos_ind)
        
        mm=m(mpos_ind);
        bb=b(mpos_ind);
        
        F(mpos_ind)=(1./sqrt(mm)).*(elliptic12i(asin(sqrt(mm).*sin(bb)),1./mm)); %cannot output complex part
        % disp('complex part may be missing');
    end
    
    mreg_ind=m<=1&m>=0 & ~phase_ind;
    
    if any(mreg_ind)
        
        mm=m(mreg_ind);
        bb=b(mreg_ind);
        
        F(mreg_ind)=elliptic12i(bb,mm);
        
    end

end
      

% E function

if nargout>1
if nargin==1  %now b= pi/2
    
    m=b;
    
    
    E=nan(size(m));
    
    if any(m<0)
        
        mm=m(m<0);
        
        [FF,EE]= ellipke(-mm./(1-mm)); %to define if your using the F output in elliptic12 or the E output
        E(m<0)=sqrt(1-mm).*EE;
    end
    
    if any(m>1)
        
        mm=m(m>1);
        
        [FF,EE]=elliptic12i(asin(sqrt(mm)),1./mm);      
        E(m>1)=((1./sqrt(mm))-sqrt(mm)).*FF+sqrt(mm).*EE;
    end
    
    
    if any(m<=1&m>=0)
        
        mm=m(m<=1&m>=0);
        [FF,E(m<=1&m>=0)]=ellipke(mm);
        
    end    

    
elseif nargin==2
   
   phase_ind = b>pi/2 | b<0;
   
    if any(phase_ind)
        
        mm = m(phase_ind);
        bb = b(phase_ind);
        
        phi = mod(bb+pi/2,pi)-pi/2;
        a = round(bb./pi);
        [F1,E1]=elliptic12x(mm);
        [FF,EE]=elliptic12x(abs(phi),mm);
        E(phase_ind) = 2*a.*E1 + sign(phi).*EE;
    end
    
   mz_ind = m==0 & ~phase_ind;
    if any(mz_ind)
        
        bb=b(mz_ind);
        
        E(mz_ind)=bb;
    end
    
    mneg_ind = m<0 & ~phase_ind;
    if any(mneg_ind)
        
        mm=m(mneg_ind);
        bb=b(mneg_ind);
       
       t=asin((sin(bb).*sqrt(1-mm))./sqrt(1-mm.*(sin(bb)).^2));
       [FF,EE]= elliptic12i(t,-mm./(1-mm)); %to define if your using the F output in elliptic12 or the E output
       E(mneg_ind)=mm.*(sin(t).*cos(t)./sqrt(1-mm.*(cos(t)).^2))+sqrt(1-mm).*EE;
    end
    
    mpos_ind = m>1 & ~phase_ind;
    
    if any(mpos_ind)
        
        mm=m(mpos_ind);
        bb=b(mpos_ind);
        
        [FF,EE]=elliptic12i(asin(sqrt(mm).*sin(bb)),1./mm); %cannot display complex part      
        E(mpos_ind)=((1./sqrt(mm))-sqrt(mm)).*FF+sqrt(mm).*EE;
        warning('elliptic123:BadComplex','Complex part may be missing');
        
    end

    mreg_ind=m<=1&m>0 & ~phase_ind;
    
    if any(mreg_ind)
        
        mm=m(mreg_ind);
        bb=b(mreg_ind);
        [FF,E(mreg_ind)]=elliptic12i(bb,mm); %note the elliptic12i function cannot evaluate at m=0 for E
        
    end
      
end
end
end

function [P]=elliptic3x(b,m,n) 

%   incomplete case
%   P=elliptic3(b,m,n) where b is the phase angle in the range of 0<b<Pi/2
%   module m where -inf<m<1 
%   m>1 can be computed provided that the input phase angle is not complex 
%   elliptic3 is not able to read complex inputs
%   parameter n<1 

%   For the complete case P=elliptic3(m,n)
%   Same conditions as above where b=Pi/2
%   simpler equations implemented (ellippi and ellippin)
%   n>1 possible but complex part of solution not shown



if nargin==2  %now b= pi/2 so its in the form of ellipticP(m,n)
    
    n=m;
    m=b;
    
    if length(m)==1, m=m(ones(size(n))); end
    if length(n)==1, n=n(ones(size(m))); end
    
    if ~isequal(size(n),size(m))
        error('m and b must be equal sizes')
    end

    mz_ind = m==0;
    
    if any(mz_ind)
        
        nn=n(mz_ind);
        
        P(mz_ind)=pi./(2.*sqrt(1-nn));
        
    end 
    
    mlnl_ind = m~=0 & m<=1 & n<=1; % mless than 1 and nless than 1 index
    
    if any(mlnl_ind)
        
        nn=n(mlnl_ind);
        mm=m(mlnl_ind);
        
        P(mlnl_ind)=ellippi(nn,mm);
        
        %P=1/((n-m)*sqrt(1-m))*(-m*ellipke(-m/(1-m))+n*ellippi((n-m)/(1-m),-m/(1-m)));
        
    end 
    
    mgnl_ind=m>1 & n<=1;
    
    if any(mgnl_ind) %% this doesnt work as b becomes complex
        
        mm=m(mgnl_ind);
        nn=n(mgnl_ind);
        
        P(mgnl_ind)=1./sqrt(mm).*elliptic3ic(asin(sqrt(mm)),1./mm,nn./mm);
     
    end
    
    mlng_ind=n>1 & m<=1;
    
    if any(mlng_ind) %%only gives real
        
        mm=m(mlng_ind);
        nn=n(mlng_ind);
        
        P(mlng_ind)=ellippin(nn,mm);
        
        warning('elliptic123:MissingComplex','Cannot compute complex part for Elliptic3 with m<=1 and n>1')
        
    end
    
    mgng_ind=n>1 & m>1;
    
    if any(mgng_ind)
            
        error('n and m vectors cannot have 1 component greater than 1')
            
    end
            
elseif nargin==3
    
    
    isize=max(max(size(b),size(m)),size(n));
    
    
    if length(n)==1, n=n(ones(isize)); end
    if length(m)==1, m=m(ones(isize)); end
    if length(b)==1, b=b(ones(isize)); end
    
    if ~isequal(size(n),size(m))|| ~isequal(size(n),size(b))
        
        error('m and b and n must be equal sizes')
    end
    
    phase_ind = b>pi/2 | b<0;
    mone_ind= m==1;
    
    phasen_ind=phase_ind & n>1;
    phasenm_ind=phasen_ind & mone_ind; 
    
    if any(phasen_ind)
        
        mm = m(phasen_ind);
        bb = b(phasen_ind);
        nn = n(phasen_ind);
   
    
        if any(b>pi/2 & b<pi) 
            
            cc=bb-pi/2;
            
            P(phasen_ind)=conj(-elliptic3x(pi/2-cc,mm,nn));
        
        elseif any(b>pi)
            
            %a=round(bb./pi);
            
            P(phasen_ind)=conj(elliptic3x(bb-pi,mm,nn)); 
        
        elseif any(b<-pi/2 & b>-pi)
            
            cc=-pi/2-bb;
            
            P(phasen_ind)=conj(-elliptic3x(-pi/2+cc,mm,nn));
            
        elseif any(b<pi)
            
            P(phasen_ind)=conj(elliptic3x(bb+pi,mm,nn));
        end
        
    end

        %if any(phasenm_ind)
            
         %P(phasenm_ind)=-elliptic3x(bb,mm,nn);

        
        %end 
         %if any(phasenm_ind)
             
             %disp('fff')
             %mm = m(phasenm_ind);
             %bb = b(phasenm_ind);
             %nn = n(phasenm_ind);
             
             %P(phasenm_ind)=(1./(1-nn)).*(log(tan(bb)+sec(bb))-0.5.*sqrt(nn).*log(1+sqrt(nn).*sin(bb))./(1-sqrt(nn).*sin(bb)));
        
         %end
         
    if any(phase_ind & ~mone_ind & n<1) %& n<1 & n>-1 When b is out side of the normal range n should be between -1 and 1 but works for n<0 anyway
    %doesnt work for n>1
        
        mm = m(phase_ind);
        bb = b(phase_ind);
        nn = n(phase_ind); 
        
        phi = mod(bb+pi/2,pi)-pi/2;
        a = round(bb./pi);
        P(phase_ind) = 2.*a.*elliptic3x(mm,nn) + sign(phi).*elliptic3x(abs(phi),mm,nn);
    end
    
    if any(phase_ind & mone_ind & ~phasen_ind)
        
        P(mone_ind)=inf;
    end 
    
    M=(1./sin(b)).^2; %critical value which goes from real inputs to complex inputs
    
    m_ind= m>M;
    
    if any(m_ind)
        
        error('one of your m value inputs is greater than the critical value')
        
    end
    
    mnequal_ind = m==n & ~phase_ind;
    
    if any(mnequal_ind)
        
        bb=b(mnequal_ind);
        mm=m(mnequal_ind);
        nn=n(mnequal_ind);
        
        [FF,EE]= elliptic12x(bb,mm);
        
        P(mnequal_ind)=(1./(1-mm)).*(EE-((mm./sqrt(1-mm.*(sin(bb)).^2)).*sin(bb).*cos(bb)));
        
    end
   
    mgnl_ind = m>1 & n<1 & ~phase_ind;
    
    if any(mgnl_ind)
        
        bb=b(mgnl_ind);
        mm=m(mgnl_ind);
        nn=n(mgnl_ind);
        
        P(mgnl_ind)=1./sqrt(mm).*elliptic3ic(asin(sqrt(mm).*sin(bb)),1./mm,nn./mm);
        
    end 
        
    mlnl_ind = m~=n & m<0 & n<1 & ~phase_ind;
    
    if any(mlnl_ind)
        
        bb=b(mlnl_ind);
        mm=m(mlnl_ind);
        nn=n(mlnl_ind);
        
        t=asin((sin(bb).*sqrt(1-mm))./sqrt(1-mm.*(sin(bb)).^2));
        
        P(mlnl_ind)=1./((nn-mm).*sqrt(1-mm)).*(-mm.*elliptic12x(t,-mm./(1-mm))+nn.*elliptic3ic(t,-mm./(1-mm),(nn-mm)./(1-mm)));
    
    end

    mnormnl_ind=n<1 & ~phase_ind & m>=0 & m<=1; %when m is between [0 1] and n<1
     
    if any(mnormnl_ind)
        
        bb=b(mnormnl_ind);
        mm=m(mnormnl_ind);
        nn=n(mnormnl_ind);
        
        %P(mnormnl_ind)=(1./(1-nn)).*(log(tan(bb)+sec(bb))-0.5.*sqrt(nn).*log(1+sqrt(nn).*sin(bb))./((1-sqrt(nn).*sin(bb))));
        P(mnormnl_ind)=elliptic3ic(bb,mm,nn);
        
    end
    
    ng_ind=n>1 & m<n & ~phase_ind; %case where n>1 but m<n
    
    if any(ng_ind)
        
        bb=b(ng_ind);
        mm=m(ng_ind);
        nn=n(ng_ind);
        
        N=mm./nn;
        P1=sqrt(((nn-1).*(1-mm./nn))); %refer to 17.7.8 in abramowitz
        D=sqrt(1-mm.*(sin(bb)).^2);
        
        %P(ng_ind)=(1./(1-nn)).*(log(tan(bb)+sec(bb))-0.5.*sqrt(nn).*log(1+sqrt(nn).*sin(bb))./((1-sqrt(nn).*sin(bb))));
        P(ng_ind)=-elliptic3x(bb,mm,N)+elliptic12x(bb,mm)+(1./(2.*P1)).*log((D+P1.*tan(bb)).*(D-P1.*tan(bb)).^-1);
       
    end
    
    ngreat_ind=n>1 & ~phase_ind & m>n;
        
    if any(ngreat_ind)
        
        error('one of your n inputs is not less than 1')
    
    end

  
end

end

function [Fi,Ei,Zi] = elliptic12i(u,m,tol)

% ELLIPTIC12i evaluates the Incomplete Elliptic Integrals 
% of the First, Second Kind and Jacobi's Zeta Function for the complex 
% value of phase U. Parameter M must be in the range 0 <= M <= 1. 
%
%   [Fi,Ei,Zi] = ELLIPTIC12i(U,M,TOL) where U is a complex phase in 
%   radians, M is the real parameter and TOL is the tolerance (optional). 
%   Default value for the tolerance is eps = 2.220e-16.
%
%   ELLIPTIC12i uses the function ELLIPTIC12 to evaluate the values of
%   corresponding integrals.
%
%   Example:
%   [phi1,phi2] = meshgrid(-2*pi:3/20:2*pi, -2*pi:3/20:2*pi);
%   phi = phi1 + phi2*i;
%   [Fi,Ei,Zi] = elliptic12i(phi, 0.5);
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC12.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html 
% Everyone is permitted to copy and distribute verbatim copies of this 
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE. 
%  
% Copyright (C) 2007 by Moiseev Igor. All rights reserved.
% 34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
% For support, please reply to 
%     moiseev.igor[at]gmail.com, moiseev[at]sissa.it
%     Moiseev Igor, 
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(m)
    error('The parameter M must be real.')
end

if any(m < 0) || any(m > 1) 
    error('M must be in the range 0 <= M <= 1.'); 
end

% if the input is real, evaluate the elliptic integrals with ELLIPTIC12
% if isreal(u)
%    [Fi,Ei,Zi] = elliptic12(u,m,tol);
%    return;
% end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u))
    error('U and M must be the same size.'); 
end

% capture memory and save the structure of input arrays
F1 = zeros(size(u)); F2 = zeros(size(u)); 
E1 = F1;     E2 = F1;
Z1 = F1;     Z2 = F1;
Fi = F1;     Ei = F1;
Zi = F1;
lambda = []; mu = []; 
I = [];      J  = [];

% make a row vector
m = m(:).'; 
u = u(:).';

% represent u in the form u = phi + i*psi
phi = real(u);
psi = imag(u);

% to avoid singularity of COT(phi) at zero add EPS
I = find (abs(phi) < eps);
phi(I) = eps;
I = [];

% finding the roots of the equation
% X^2 - (cot(phi)^2+m*sinh(psi)^2*csc(phi)^2-1+m)X - (1-m)*cot(phi)^2 = 0
b = -(cot(phi).^2 + m.*sinh(psi).^2.*csc(phi).^2-1+m);
c = -(1-m).*cot(phi).^2;

X1 = -b/2 + sqrt(b.^2/4-c);
I = find(X1>=0);

if length(I) ~= length(u)
    X2 = -b/2 - sqrt(b.^2/4-c);
    J = find(X2>=0);
end

if( ~isempty(I) ) 
    lambda(I) = acot( sqrt(X1(I)) ); 
    mu(I)     = atan( sqrt(1./m(I).*(tan(phi(I)).^2.*cot(lambda(I)).^2 - 1)) );
end
if( ~isempty(J) ) 
    lambda(J) = acot( sqrt(X2(J)) ); 
    mu(J)     = atan( sqrt(1./m(J).*(tan(phi(J)).^2.*cot(lambda(J)).^2 - 1)) );
end

% change of variables taking into account periodicity ceil to the right
lambda = (-1).^floor(phi/pi*2).*lambda + pi*ceil(phi/pi-0.5+eps);
mu     = sign(psi).*real(mu);

[F1(:),E1(:)] = elliptic12ic(lambda, m, tol);
[F2(:),E2(:)] = elliptic12ic(mu, 1-m, tol);
 
% complex values of elliptic integral of the first kind
Fi = F1 + sqrt(-1)*F2;

% some calucation optimiziation
sin_lam = sin(lambda); cos_lam = cos(lambda);
sin_mu = sin(mu); cos_mu = cos(mu);

b1 = m.*sin_lam.*cos_lam.*sin_mu.^2.*sqrt(1-m.*sin_lam.^2);
b2 = sin_mu.*cos_mu.*(1-m.*sin_lam.^2).*sqrt(1-(1-m).*sin_mu.^2);
b3 = cos_mu.^2 + m.*sin_lam.^2.*sin_mu.^2;

% complex values of elliptic integral of the second kind
Ei(:) = (b1 + sqrt(-1)*b2)./b3;
Ei(:) = Ei(:) + E1(:) + sqrt(-1)*(-E2(:) + F2(:));

[K,Ee] = ellipke(m);
% complex values of zeta function
Zi(:) = Ei(:) - Ee(:)./K(:).*Fi(:);
end

% END FUNCTION ELLIPTIC12i()

function [F,E,Z] = elliptic12ic(u,m,tol)

% ELLIPTIC12ic evaluates the value of the Incomplete Elliptic Integrals 
% of the First, Second Kind and Jacobi's Zeta Function.
%
%   [F,E,Z] = ELLIPTIC12(U,M,TOL) where U is a phase in radians, 0<M<1 is 
%   the module and TOL is the tolerance (optional). Default value for 
%   the tolerance is eps = 2.220e-16.
%
%   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean 
%   and Descending Landen Transformation described in [1] Ch. 17.6,
%   to determine the value of the Incomplete Elliptic Integrals 
%   of the First, Second Kind and Jacobi's Zeta Function [1], [2].
%
%       F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
%       E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
%       Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
%
%   Tables generating code ([1], pp. 613-621):
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees
%       [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC12I, ELLIPTIC3, THETA, AGM.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html 
% Everyone is permitted to copy and distribute verbatim copies of this 
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE. 
%  
% Copyright (C) 2007 by Moiseev Igor. All rights reserved.
% 34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
% For support, please reply to 
%     moiseev[at]sissa.it, moiseev.igor[at]gmail.com
%     Moiseev Igor, 
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
%
% The code is optimized for ordered inputs produced by the functions 
% meshgrid, ndgrid. To obtain maximum performace (up to 30%) for singleton, 
% 1-dimensional and random arrays remark call of the function unique(.) 
% and edit further code. 

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) || ~isreal(m)
    error('Input arguments must be real. Use ELLIPTIC12i for complex arguments.');
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

F = zeros(size(u)); 
E = F;              
Z = E;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

I = uint32( find(m ~= 1 & m ~= 0) );
if ~isempty(I)
    [mu,J,K] = unique(m(I));   % extracts unique values from m
    K = uint32(K);
    mumax = length(mu);
    signU = sign(u(I));

    % pre-allocate space and augment if needed
        chunk = 7;
        a = zeros(chunk,mumax);
        c = a; 
        b = a;
        a(1,:) = ones(1,mumax);
        c(1,:) = sqrt(mu);
        b(1,:) = sqrt(1-mu);
        n = uint32( zeros(1,mumax) );
        i = 1;
        while any(abs(c(i,:)) > tol)                                    % Arithmetic-Geometric Mean of A, B and C
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,mumax)];
          b = [b; zeros(2,mumax)];
          c = [c; zeros(2,mumax)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        in = uint32( find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol)) );
        if ~isempty(in)
          [mi,ni] = size(in);
          n(in) = ones(mi,ni)*(i-1);
        end
        end
     
    mmax = length(I);
        mn = double(max(n));
        phin = zeros(1,mmax);     C  = zeros(1,mmax);    
        Cp = C;  e  = uint32(C);  phin(:) = signU.*u(I);
        i = 0;   c2 = c.^2;
        while i < mn                                                    % Descending Landen Transformation 
        i = i + 1;
        in = uint32(find(n(K) > i));
        if ~isempty(in)     
            phin(in) = atan(b(i,K(in))./a(i,K(in)).*tan(phin(in))) + ...
                pi.*ceil(phin(in)/pi - 0.5) + phin(in);
            e(in) = 2.^(i-1) ;
            C(in) = C(in)  + double(e(in(1)))*c2(i,K(in));
            Cp(in)= Cp(in) + c(i+1,K(in)).*sin(phin(in));  
        end
        end
    
    Ff = phin ./ (a(mn,K).*double(e)*2);                                                      
    F(I) = Ff.*signU;                                               % Incomplete Ell. Int. of the First Kind
    Z(I) = Cp.*signU;                                               % Jacobi Zeta Function
    E(I) = (Cp + (1 - 1/2*C) .* Ff).*signU;                         % Incomplete Ell. Int. of the Second Kind
end

% Special cases: m == {0, 1}
m0 = find(m == 0);
if ~isempty(m0), F(m0) = u(m0); E(m0) = u(m0); Z(m0) = 0; end

m1 = find(m == 1);
um1 = abs(u(m1)); 
if ~isempty(m1), 
    N = floor( (um1+pi/2)/pi );  
    M = find(um1 < pi/2);              
    
    F(m1(M)) = log(tan(pi/4 + u(m1(M))/2));   
    F(m1(um1 >= pi/2)) = Inf.*sign(u(m1(um1 >= pi/2)));
    
    E(m1) = ((-1).^N .* sin(um1) + 2*N).*sign(u(m1)); 
    
    Z(m1) = (-1).^N .* sin(u(m1));                      
end
end 


function Pi = elliptic3ic(u,m,c)

if nargin<3, error('Not enough input arguments.'); end
if ~isreal(u) | ~isreal(m) | ~isreal(c)
    error('Input arguments must be real.')
end
if any(m < 0) | any(m > 1),  
  error('M must be in the range [0, 1].');
end
if any(c > 1),  
  error('C must be in the range [-inf, 1].');
end
if any(u > pi/2) | any(u < 0),  
    error('U must be in the range [0, pi/2].'); 
end

[mm,nm] = size(m);
[mu,nu] = size(u);
if length(m)==1, m = m(ones(size(u))); end
if length(c)==1, c = c(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m), size(c), size(u)), 
        error('U, M and C must be the same size.'); 
end

Pi = zeros(size(u));
m = m(:).';    % make a row vector
u = u(:).';
c = c(:).';

I = find( u==pi/2 & m==1 | u==pi/2 & c==1 );

t = [ 0.9931285991850949,  0.9639719272779138,...            % Base points 
      0.9122344282513259,  0.8391169718222188,...            % for Gauss-Legendre integration
      0.7463319064601508,  0.6360536807265150,...
      0.5108670019508271,  0.3737060887154195,...
      0.2277858511416451,  0.07652652113349734 ];                             
w = [ 0.01761400713915212, 0.04060142980038694,...           % Weights
      0.06267204833410907, 0.08327674157670475,...           % for Gauss-Legendre integration
      0.1019301198172404,  0.1181945319615184,...
      0.1316886384491766,  0.1420961093183820,...
      0.1491729864726037,  0.1527533871307258  ];
  
P = 0;  i = 0;
while i < 10
    i  = i + 1;
    c0 = u.*t(i)/2;
    P  = P + w(i).*(g(u/2+c0,m,c) + g(u/2-c0,m,c));
end
P = u/2.*P;
Pi(:) = P;                                                   % Incomplete elliptic integral of the third kind

% special values u==pi/2 & m==1 | u==pi/2 & c==1
Pi(I) = inf;
return;


function g = g(u,m,c)
%  g = 1/((1 - c*sin(u)^2)*sqrt(1 - m*sin(u)^2));

 sn2 = sin(u).^2;
 g = 1./((1 - c.*sn2).*sqrt(1 - m.*sn2));
return;

end

end

function PI = ellippi(n,m)

% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: http://dlmf.nist.gov/19.8
%
% Valid for n <= 1 and m <= 1


a0 = 1;
g0 = sqrt(1-m);
s0 = m;
nn = 0;

p0 = sqrt(1-n);
Q0 = 1;
QQ = Q0;

w1 = ones(size(m));

while max(w1(:)) > eps % assume ellip2 converges slower than ellip3

 % for Elliptic I
 a1 = (a0+g0)/2;
 g1 = sqrt(a0.*g0);

 % for Elliptic II
 nn = nn + 1;
 c1 = (a0-g0)/2;
 w1 = 2^nn*c1.^2;
 s0 = s0 + w1;

 % for Elliptic III
 rr = p0.^2+a0.*g0;
 p1 = rr./(2.*p0);
 Q1 = 0.5*Q0.*(p0.^2-a0.*g0)./rr;
 QQ = QQ+Q1;

 a0 = a1;
 g0 = g1;
 Q0 = Q1;
 p0 = p1;

end

PI = pi./(4.*a1).*(2+n./(1-n).*QQ);

im = find(m == 1);
if ~isempty(im)
 k(im) = inf;
 e(im) = ones(length(im),1);
 PI(im) = inf;
end

end

function PI = ellippin(n,m)

% Complete elliptic integrals calculated with the arithmetric-geometric mean
% algorithms contained here: http://dlmf.nist.gov/19.8
%
% Valid for n > 1 and m <= 1


a0 = 1;
g0 = sqrt(1-m);
s0 = m;
nn = 0;

p0 = sqrt(1-(m/n));
Q0 = 1;
QQ = Q0;

w1 = ones(size(m));

while max(w1(:)) > eps % assume ellip2 converges slower than ellip3

 % for Elliptic I
 a1 = (a0+g0)/2;
 g1 = sqrt(a0.*g0);

 % for Elliptic II
 nn = nn + 1;
 c1 = (a0-g0)/2;
 w1 = 2^nn*c1.^2;
 s0 = s0 + w1;

 % for Elliptic III
 rr = p0.^2+a0.*g0;
 p1 = rr./(2.*p0);
 Q1 = 0.5*Q0.*(p0.^2-a0.*g0)./rr;
 QQ = QQ+Q1;

 a0 = a1;
 g0 = g1;
 Q0 = Q1;
 p0 = p1;

end

PI = pi./(4.*a1).*((m/(m-n)).*QQ);

im = find(m == 1);
if ~isempty(im)
 k(im) = inf;
 e(im) = ones(length(im),1);
 PI(im) = inf;
end

end

