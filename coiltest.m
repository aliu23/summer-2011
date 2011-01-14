displ=linspace(0,0.02,100);
f=displ;
for ii=1:length(displ)
    f(ii)=magnetcoil(displ(ii),9e-3,10e-3,1,10e-3,10.5e-3,20e-3,100,1);
end
plot(displ,f)
