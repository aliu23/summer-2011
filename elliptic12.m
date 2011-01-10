% Elliptic 12 computes the First and Second Elliptic integrals for both the
% complete and incomplete cases. 

%if in the form [F,E]=elliptic12(b,m)
%outputs F and E integrals with phase angle b where 0<b<pi/2

%For the incomplete Case: 
%when complex numbers are expected, the complex part will not be outputed

%For complete case the form is [F,E]=elliptic12(m) where b is Pi/2
%m can take all values here
%will output the whole answer



function [F,E]=elliptic12(b,m)

%F function

 
if nargin==1  %now b= pi/2
    
    m=b;
    
    if m<0 
    
        
        F=(1/sqrt(1-m))*ellipke(-m/(1-m)); %in future may want to consider complex inputs
    
    elseif m>1
        
        F=(1/sqrt(m))*(elliptic12i(asin(sqrt(m)),1/m)); 
        

    else %0<m<=1
          
        F=ellipke(m); 
        
    end
    
elseif nargin==2
    
    if b>pi/2 || b<0
        
        phi = mod(b+pi/2,pi)-pi/2;
        a = round(b/pi);
        F = 2*a.*elliptic12(m) + sign(phi)*elliptic12(abs(phi),m);
   
    elseif m<0 
    
        t=asin((sin(b)*sqrt(1-m))/sqrt(1-m*(sin(b))^2)); %t is actually theta
        F=(1/sqrt(1-m))*elliptic12i(t,-m/(1-m)); 
    
    
    elseif m>1 
        
        F=(1/sqrt(m))*(elliptic12i(asin(sqrt(m)*sin(b)),1/m)); %cannot output complex part
        disp('complex part may be missing');
        
    else  
          
        F=elliptic12i(b,m);
        
    
    end
      
end


% E function


if nargin==1  %now b= pi/2
    
    m=b;
    
    if m<0
           
        [FF,EE]= ellipke(-m/(1-m)); %to define if your using the F output in elliptic12 or the E output
        E=sqrt(1-m)*EE;
    
    elseif m>1
        
        [FF,EE]=elliptic12i(asin(sqrt(m)),1/m);      
        E=((1/sqrt(m))-sqrt(m))*FF+sqrt(m)*EE;
     
    else
        
        
        [FF,E]=ellipke(m);
        
        
    end
    
elseif nargin==2
   
    if b>pi/2 || b<0
        
        phi = mod(b+pi/2,pi)-pi/2;
        a = round(b/pi);
        [F1,E1]=elliptic12(m);
        [FF,EE]=elliptic12(abs(phi),m);
        E = 2*a.*E1 + sign(phi)*EE;
        
   elseif m<0
        t=asin((sin(b)*sqrt(1-m))/sqrt(1-m*(sin(b))^2));   
        [FF,EE]= elliptic12(t,-m/(1-m)); %to define if your using the F output in elliptic12 or the E output
        E=m*(sin(t)*cos(t)/sqrt(1-m*(cos(t))^2))+sqrt(1-m)*EE;
    
 elseif m>1 
        
        [FF,EE]=elliptic12i(asin(sqrt(m)*sin(b)),1/m); %cannot display complex part      
        E=((1/sqrt(m))-sqrt(m))*FF+sqrt(m)*EE;
        disp('complex part may be missing');
     
 else 
        
        
        [FF,E]=elliptic12i(b,m);
        
        
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
