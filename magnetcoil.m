function [ F ] = magnetcoil(z,rm,lm,Br,rc,Rc,lc,N,I)
%MAGNETCOIL Function which evaluates the force between a current carrying
%coil and a magnet
%
% Both coil and magnet are cylindrical with the following dimensions:
%
%   z = axial displacement between coil/magnet centres (m)
%
%  rm = magnet radius (m)
%  lm = magnet length (m)
%  Br = magnetisation strength (T)
%
%  rc = inner coil radius (m)
%  Rc = outer coil radius (m)
%  lc = coil length (m)
%   N = number of coil turns
%   I = coil current (A)
%
% The force is that acted upon the magnet and the coil is assumed to be
% wound in the anti-clockwise direction when viewed from the end in the
% positive axial-direction. (I.e., the coil will attact the magnet.)

Q = nan(size(z));

for ii = 1:length(z)
  Q(ii)=dblquad(@(rr,zz) auxcoil(z(ii),rr,zz),...
    rc,Rc,-lc./2,lc./2);
end
  
F=(N.*I.*Br./(lc.*(Rc-rc))).*Q;


  function [F]=auxcoil(z,r1,z1)
    
    m3=2.*rm.*r1;
    
    m2p=sqrt((rm+r1).^2+(z+lm./2-z1).^2);
    m2n=sqrt((rm+r1).^2+(z-lm./2-z1).^2);
    
    m1p=(2.*m3)./(m2p.^2);
    m1n=(2.*m3)./(m2n.^2);
    
    [Fp,Ep]=ellipke(m1p);
    [Fn,En]=ellipke(m1n);
    
    fzp =  (m2p-(m3./m2p)).*Fp-m2p.*Ep;
    fzn = -(m2n-(m3./m2n)).*Fn+m2n.*En;
    
    F=fzp+fzn;
    
  end

end
