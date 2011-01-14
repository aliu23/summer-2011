
%MAGNETCOIL Function which evaluates the force between a current carrying
%coil and a magnet
%z=magnet length
%rm=magnet radius
%lm=magnet length
%Br=magnetisation
%rc=inner coil radius
%Rc=outer coil radius
%lc=coil length
%N=coil turns
%I=coil current

function [ F ] = magnetcoil(z,rm,lm,Br,rc,Rc,lc,N,I)

Q=dblquad(@auxcoil,rc,Rc,-lc./2,lc./2); %evaluate the double integral

F=(N.*I.*Br./(lc.*(Rc-rc))).*Q;



    function [F]=auxcoil(r1,z1)
        

zetap=z+lm./2; %when ei is positve

zetan=z-lm./2; %when ei is negative

m3=2.*rm.*r1;

m2p=sqrt((rm+r1).^2+(zetap-z1).^2);

m2n=sqrt((rm+r1).^2+(zetan-z1).^2);

m1p=(2.*m3)./(m2p.^2);

m1n=(2.*m3)./(m2n.^2);

[Fp,Ep]=elliptic12(m1p);

[Fn,En]=elliptic12(m1n);

fzp = (m2p-(m3./m2p)).*Fp-m2p.*Ep; %positive case of fz

fzn = -((m2n-(m3./m2n)).*Fn-m2n.*En);

F=fzp+fzn;

    end

end


