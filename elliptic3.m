function [P]=elliptic3(b,m,n) 

% ELLIPTIC3 evaluates incomplete elliptic integral of the third kind.
%   Pi = ELLIPTIC3(U,M,C) where U is a phase in radians, 0<M<1 is 
%   the module and 0<C<1 is a parameter. 
%
%   ELLIPTIC3 uses Gauss-Legendre 10 points quadrature template 
%   described in [3] to determine the value of the Incomplete Elliptic 
%   Integral of the Third Kind (see [1, 2]).
%
%   Pi(u,m,c) = int(1/((1 - c*sin(t)^2)*sqrt(1 - m*sin(t)^2)), t=0..u)
%
%   Tables generating code ([1], pp. 625-626):
%	    [phi,alpha,c] = meshgrid(0:15:90, 0:15:90, 0:0.1:1); 
%   	Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
%  
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 17.7.
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989.
%   [3] S. Zhang, J. Jin "Computation of Special Functions" (Wiley, 1996).

%   For support, please reply to 
%       moiseev[at]sissa.it
%       Moiseev Igor, 
%       34106, SISSA, via Beirut n. 2-4,  Trieste, Italy

%   additions to script added 6/01/2011

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
        
        %disp('complex part may be missing')
        
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
    
    if any(phase_ind & ~mone_ind) %& n<1 & n>-1 When b is out side of the normal range n should be between -1 and 1 but works for n<0 anyway
    %doesnt work for n>1
        
        mm = m(phase_ind);
        bb = b(phase_ind);
        nn = n(phase_ind); 
        
        phi = mod(bb+pi/2,pi)-pi/2;
        a = round(bb./pi);
        P(phase_ind) = 2.*a.*elliptic3(mm,nn) + sign(phi).*elliptic3(abs(phi),mm,nn);
    end
    
    if any(phase_ind & mone_ind)
        
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
        
        [FF,EE]= elliptic12(bb,mm);
        
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
        
        P(mlnl_ind)=1./((nn-mm).*sqrt(1-mm)).*(-mm.*elliptic12(t,-mm./(1-mm))+nn.*elliptic3ic(t,-mm./(1-mm),(nn-mm)./(1-mm)));
    
    end

    mnormnl_ind=n<1 & ~phase_ind & m>=0 & m<=1; %when m is between [0 1] and n<0
     
    if any(mnormnl_ind)
        
        bb=b(mnormnl_ind);
        mm=m(mnormnl_ind);
        nn=n(mnormnl_ind);
        
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
        
        P(ng_ind)=-elliptic3(bb,mm,N)+elliptic12(bb,mm)+(1./(2.*P1)).*log((D+P1.*tan(bb)).*(D-P1.*tan(bb)).^-1);
       
    end
    
    ngreat_ind=n>1 & ~phase_ind & m>n;
        
    if any(ngreat_ind)
        
        error('one of your n inputs is not less than 1')
    
    end
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
