%% Test suite for additions to elliptic12 and elliptic3

clc

m = -5:0.5:5;

%% Elliptic F and E

fprintf('Tests for EllipticF and EllipticE\n=================================\n\n')

%% Complete F and E

ellKvalues = [0.9555, 0.981, 1.0095, 1.0416, 1.0783, 1.1209, 1.1714, 1.233, 1.311, 1.4157, 1.5708, 1.8541, inf, 1.6566 - 1.4157*i, 1.311 - 1.311*i, 1.1242 - 1.233*i, 1.0011 - 1.1714*i, 0.9117 - 1.1209*i, 0.8429 - 1.0783*i, 0.7877 - 1.0416*i, 0.7422 - 1.0095*i];
ellEvalues = [2.8302, 2.7347, 2.6352, 2.5312, 2.4221, 2.3069, 2.1844, 2.053, 1.9101, 1.7518, 1.5708, 1.3506, 1., 0.7163 + 0.336*i, 0.5991 + 0.5991*i, 0.5263 + 0.82*i, 0.4752 + 1.013*i, 0.4367 + 1.186*i, 0.4063 + 1.3439*i, 0.3815 + 1.4897*i, 0.3608 + 1.6257*i];

[tryk tryE] = elliptic12(m);
tryk = round(10000*tryk)/10000;
tryE = round(10000*tryE)/10000;

fprintf('Complete F for all m: ')
if all(tryk==ellKvalues)
  fprintf('passed\n')
else
  fprintf('failed\n')
end

fprintf('Complete E for all m: ')
if all(tryE==ellEvalues) 
  fprintf('passed\n')
else
  fprintf('failed\n')
end


%% Incomplete b > pi/2

fprintf('\nIncomplete tests\n----------------\n\nPhase between 0 and pi/2:\n\n')

b=0.5;
M = (1./sin(b)).^2;
disp(['The critical value of m for b=',num2str(b),' is M=',num2str(M)])

%% incomplete F elliptic12i cannot output the imaginary part for m values
% M and over,, real part ok
 
ellFvalues1 = [0.4325, 0.4374, 0.4426, 0.4481, 0.454, 0.4602, 0.4669, 0.4742, 0.482, 0.4906, 0.5, 0.5105, 0.5222, 0.5357, 0.5514, 0.5702, 0.5938, 0.6258, 0.6774, 0.7877 - 0.0986*i, 0.7422 - 0.1898*i];

tryF = elliptic12(b,m);
tryF = round(10000*tryF)/10000;

fprintf('incomplete F(m<=1,0<b<pi/2): ')

if all(tryF(m<=1)==ellFvalues1(m<=1))
  fprintf('passed\n')
elseif all(real(tryF(m<=1))==real(ellFvalues1(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete F(m<M, 0<b<pi/2): ')

if all(tryF(m<M)==ellFvalues1(m<M)) 
  fprintf('passed\n')
elseif all(real(tryF(m<M))==real(ellFvalues1(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')  
end

fprintf('Incomplete F(m>M, 0<b<pi/2): ')

if all(tryF(m>M)==ellFvalues1(m>M)) 
  fprintf('passed\n')
elseif all(real(tryF(m>M))==real(ellFvalues1(m>M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

%% incomplete E (same prob as F)

ellEvalues1 = [0.5864, 0.5787, 0.5707, 0.5627, 0.5544, 0.5459, 0.5372, 0.5283, 0.5192, 0.5097, 0.5, 0.4899, 0.4794, 0.4685, 0.457, 0.4448, 0.4319, 0.4177, 0.4018, 0.3815 + 0.0011*i, 0.3608 + 0.0093*i];

[F,tryE]=elliptic12(b,m);
tryE = round(10000*tryE)/10000;

fprintf('Incomplete E(m<=1,0<b<pi/2): ')

if all(tryE(m<=1)==ellEvalues1(m<=1)) 
  fprintf('passed\n')
elseif all(tryE(m<=1)==ellEvalues1(m<=1)) || all(tryE(m==0)~=ellEvalues1(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryE(m<=1))==real(ellEvalues1(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete E(m<M, 0<b<pi/2): ')

if all(tryE(m<M)==ellEvalues1(m<M)) 
  fprintf('passed\n')
  
elseif all(tryE(m<M)==ellEvalues1(m<M)) || all(tryE(m==0)~=ellEvalues1(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryE(m<M))==real(ellEvalues1(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete E(m>M, 0<b<pi/2): ')

if all(tryE(m>M)==ellEvalues1(m>M))
  fprintf('passed\n')
elseif all(real(tryE(m>M))==real(ellEvalues1(m>M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

%% b > pi/2

fprintf('\nPhase greater than pi/2:\n\n')

b=2;
M = (1./sin(b)).^2;
disp(['The critical value of m for b=',num2str(b),' is M=',num2str(M)])

%% incomplete F b>pi/2 (for all m<M works) 
%for m>M correct up until we encounter complex outputs then it just spits
%out the reals

ellFvalues2 = [1.1354, 1.1688, 1.2063, 1.2489, 1.2979, 1.3554, 1.4244, 1.5095, 1.6192, 1.7697, 2., 2.4444, inf, 1.6566 - 2.0956*i, 1.311 - 1.7707*i, 1.1242 - 1.6035*i, 1.0011 - 1.4903*i, 0.9117 - 1.405*i, 0.8429 - 1.337*i, 0.7877 - 1.2807*i, 0.7422 - 1.2329*i];
tryF = nan(size(m));

tryF = elliptic12(b,m);
tryF = round(10000*tryF)/10000;

fprintf('incomplete F(m<=1,b>pi/2): ')

if all(tryF(m<=1)==ellFvalues2(m<=1))
  fprintf('passed\n')
elseif all(real(tryF(m<=1))==real(ellFvalues2(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete F(m<M, b>pi/2): ')

if all(tryF(m<M)==ellFvalues2(m<M)) 
  fprintf('passed\n')
elseif all(real(tryF(m<M)==real(ellFvalues2(m<M))))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete F(m>M, b>pi/2): ')

if all(tryF(m>M)==ellFvalues2(m>M)) 
  fprintf('passed\n')
elseif all(real(tryF(m>M))==real(ellFvalues2(m>M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

%% incomplete E b>pi/2 

ellEvalues2 = [3.855, 3.7163, 3.5717, 3.4203, 3.2611, 3.0926, 2.9129, 2.7194, 2.508, 2.2722, 2., 1.6629, 1.0907, 0.7163 + 0.6099*i, 0.5991 + 1.0013*i, 0.5263 + 1.3184*i, 0.4752 + 1.5919*i, 0.4367 + 1.8354*i, 0.4063 + 2.0568*i, 0.3815 + 2.261*i, 0.3608 + 2.4513*i];

[F,tryE]=elliptic12(b,m);
tryE = round(10000*tryE)/10000;

fprintf('incomplete E(m<=1,b>pi/2): ')

if all(tryE(m<=1)==ellEvalues2(m<=1))
  fprintf('passed\n')
elseif all(real(tryE(m<=1))==real(ellEvalues2(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete E(m<M, b>pi/2): ')

if all(tryE(m<M)==ellEvalues2(m<M)) 
  fprintf('passed\n')
elseif all(real(tryE(m<M)==real(ellEvalues2(m<M))))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete E(m>M, b>pi/2): ')

if all(tryE(m>M)==ellEvalues2(m>M)) 
  fprintf('passed\n')
elseif all(real(tryE(m>M))==real(ellEvalues2(m>M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end


%% Elliptic PI

fprintf('\n\nTests for EllipticPI\n====================\n')

%% Complete

fprintf('\nComplete PI\n-----------\n\n')
fprintf('Critical value where b input goes from real to complex is M=1\n')

%% complete pi ,0<n<1, input must be real, elliptic3 cannot have complex inputs
%we need elliptic3i prob at m=0

n = 0.5;

ellPIvalues1 = [1.2566, 1.2943, 1.3367, 1.3847, 1.44, 1.5048, 1.5823, 1.6779, 1.8005, 1.9679, 2.2214, 2.7013, inf, 2.0855 - 2.4344*i, 1.5327 - 2.1041*i, 1.2659 - 1.8989*i, 1.1017 - 1.754*i, 0.9879 - 1.6439*i, 0.9032 - 1.5561*i, 0.837 - 1.4838*i, 0.7835 - 1.4227*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try 
        tryPI(ii) = elliptic3(m(ii),n);
    end
end
   
tryPI = round(10000*tryPI)/10000;




fprintf('Complete PI(m<=1,0<n<1): ')

if all(tryPI(m<=1)==ellPIvalues1(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues1(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Complete PI(m>1, 0<n<1): ')

if all(tryPI(m>1)==ellPIvalues1(m>1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>1))==real(ellPIvalues1(m>1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to a complex b input into elliptic3ic \n')
end


%% Complete Pi n>1 2 problems no elliptic3i and can only output real values
%for real inputs

n = 5;

ellPIvalues2 = [0.2271 - 0.5554*i, 0.2204 - 0.5698*i, 0.2124 - 0.5854*i, 0.2028 - 0.6024*i, 0.1911 - 0.6209*i, 0.1766 - 0.6413*i, 0.1585 - 0.6638*i, 0.1354 - 0.6888*i, 0.1048 - 0.717*i, 0.0626 - 0.7488*i, -0.7854*i, -0.1092 - 0.8279*i, inf, -0.2238 - 0.469*i, -0.1689 - 0.4558*i, -0.1417 - 0.4448*i, -0.1245 - 0.4353*i, -0.1124 - 0.4269*i, -0.1033 - 0.4194*i, -0.0961 - 0.4126*i, -0.0902 - 0.4064*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
        tryPI(ii) = elliptic3(m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Complete PI(m<=1,n>1): ')

if all(tryPI(m<=1)==ellPIvalues2(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues2(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Complete PI(m>1, n>1): ')

if all(tryPI(m>1)==ellPIvalues2(m>1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>1))==real(ellPIvalues2(m>1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to a complex b input into elliptic3ic\n')
end

   
%% complete Pi n<0 (no elliptic3i prob)

n = -5;

ellPIvalues3 = [0.4717, 0.48, 0.489, 0.499, 0.5102, 0.523, 0.5377, 0.555, 0.5762, 0.6033, 0.6413, 0.7049, inf, 0.7252 - 0.2799*i, 0.6655 - 0.2928*i, 0.6275 - 0.3023*i, 0.5984 - 0.3094*i, 0.5744 - 0.3148*i, 0.5539 - 0.3189*i, 0.5359 - 0.322*i, 0.5198 - 0.3244*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Complete PI(m<=1,n<0): ')

if all(tryPI(m<=1)==ellPIvalues3(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues2(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Complete PI(m>1, n<0): ')

if all(tryPI(m>1)==ellPIvalues3(m>1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>1))==real(ellPIvalues3(m>1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to a complex b input into elliptic3ic\n')
end


%% Incomplete PI b < 1


fprintf('\nIncomplete PI\n-------------\n\n')
fprintf('Phase between 0 and pi/2:\n\n')


b=0.5;
M = (1./sin(b)).^2;
disp(['The critical value of m for b=',num2str(b),' is M=',num2str(M)])

%% incomplete Pi 0<n<1 (elliptic 3 cannot imput complex no.)

n = 0.5;

ellPIvalues11 = [0.449, 0.4543, 0.4598, 0.4657, 0.4719, 0.4786, 0.4858, 0.4936, 0.502, 0.5112, 0.5213, 0.5326, 0.5453, 0.5598, 0.5768, 0.5973, 0.623, 0.658, 0.7149, 0.837 - 0.1111*i, 0.7835 - 0.2121*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Incomplete PI(m<=1,0<b<pi/2,0<n<1): ')

if all(tryPI(m<=1)==ellPIvalues11(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues11(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m<M, 0<b<pi/2,0<n<1): ')

if all(tryPI(m<M)==ellPIvalues11(m<M)) 
  fprintf('passed\n')
elseif all(tryPI(m<M)==ellPIvalues11(m<M)) && all(tryPI(m==0)~=ellPIvalues11(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryPI(m<M))==real(ellPIvalues11(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m>M, 0<b<pi/2,0<n<1): ')

if all(tryPI(m>M)==ellPIvalues11(m>M)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>4))==real(ellPIvalues11(m>4)))
  fprintf('real components correct, missing complex\n')
elseif any(m>M)
  fprintf('failed due to a complex b input into elliptic3ic\n')
else
  fprintf('failed')
end


%% incomplete Pi when n>1 Cant do any as n>1

n = 5;

ellPIvalues12 = [0.6705 - 0.5554*i, 0.679 - 0.5698*i, 0.688 - 0.5854*i, 0.6974 - 0.6024*i, 0.7073 - 0.6209*i, 0.7178 - 0.6413*i, 0.7288 - 0.6638*i, 0.7405 - 0.6888*i, 0.7528 - 0.717*i, 0.7658 - 0.7488*i, 0.7795 - 0.7854*i, 0.7937 - 0.8279*i, 0.8084 + 0.8781*i, 0.823 + 2.8162*i, 0.8365 + 3.0418*i, 0.8466 + 3.3322*i, 0.847 + 3.7255*i, 0.8188 + 4.3018*i, 0.6791 + 5.2686*i, -0.0961 - 1.6804*i, -0.0902 - 1.3637*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Incomplete PI(m<=1,0<b<pi/2,n>1): ')

if all(tryPI(m<=1)==ellPIvalues12(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues12(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to n>1\n')
end

fprintf('Incomplete PI(m<M, 0<b<pi/2,n>1): ')

if all(tryPI(m<M)==ellPIvalues12(m<M)) 
  fprintf('passed\n')
  
elseif all(tryPI(m<M)==ellPIvalues12(m<M)) && all(tryPI(m==0)~=ellPIvalues12(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryPI(m<M))==real(ellPIvalues12(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to n>1\n')
end

fprintf('Incomplete PI(m>M, 0<b<pi/2,n>1): ')

if all(tryPI(m>M)==ellPIvalues12(m>M)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>M))==real(ellPIvalues12(m>M)))
  fprintf('real components correct, missing complex\n')
else any(m>M)&& n>1;
  fprintf('failed due to a complex b input into elliptic3ic aswell as n>1\n')
end


%% incomplete Pi n<0, 0<b<pi/2

n=-5;
b=0.5;

ellPIvalues13 = [0.3369, 0.3401, 0.3434, 0.3469, 0.3506, 0.3546, 0.3588, 0.3634, 0.3682, 0.3735, 0.3793, 0.3856, 0.3927, 0.4007, 0.4099, 0.4207, 0.4341, 0.4518, 0.4793, 0.5359 - 0.0464*i, 0.5198 - 0.0927*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Incomplete PI(m<=1,0<b<pi/2,n<0): ')

if all(tryPI(m<=1)==ellPIvalues13(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues13(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m<M, 0<b<pi/2,n<0): ')

if all(tryPI(m<M)==ellPIvalues13(m<M)) 
  fprintf('passed\n')
elseif all(tryPI(m<M)==ellPIvalues13(m<M)) && all(tryPI(m==0)~=ellPIvalues13(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryPI(m<M))==real(ellPIvalues13(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m>M, 0<b<pi/2,n<0): ')

if all(tryPI(m>M)==ellPIvalues13(m>M)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>M))==real(ellPIvalues13(m>M)))
  fprintf('real components correct, missing complex\n')
elseif any(m>M)
  fprintf('failed due to a complex b input into elliptic3ic\n')
else
  fprintf('failed')
end 


%% Incomplete PI b > pi/2

fprintf('\nPhase greater than pi/2:\n\n')

b=2;
M = (1./sin(b)).^2;
disp(['The critical value of m for b=',num2str(b),' is M=',num2str(M)])

%% incomplete Pi 0<n<1, b>pi/2 Doesnt work for all m same no elliptic 3i problem

n = 0.5;

ellPIvalues21 = [1.5967, 1.6493, 1.7087, 1.7767, 1.8554, 1.9483, 2.0607, 2.2008, 2.3834, 2.6376, 3.0338, 3.8199, inf, 2.0855 - 3.7147*i, 1.5327 - 2.9716*i, 1.2659 - 2.5984*i, 1.1017 - 2.3562*i, 0.9879 - 2.1806*i, 0.9032 - 2.045*i, 0.837 - 1.9357*i, 0.7835 - 1.8449*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Incomplete PI(m<=1,0<b<pi/2,0<n<1): ')

if all(tryPI(m<=1)==ellPIvalues21(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues21(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m<M, 0<b<pi/2,0<n<1): ')

if all(tryPI(m<M)==ellPIvalues21(m<M)) 
  fprintf('passed\n')
  
elseif all(tryPI(m<M)==ellPIvalues21(m<M)) && all(tryPI(m==0)~=ellPIvalues21(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryPI(m<M))==real(ellPIvalues21(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m>M, 0<b<pi/2,0<n<1): ')

if all(tryPI(m>M)==ellPIvalues21(m>M)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>M))==real(ellPIvalues21(m>M)))
  fprintf('real components correct, missing complex\n')
else any(m>M);
  fprintf('failed due to a complex b input into elliptic3ic\n')
end

%% incomplete Pi n>1, b>pi/2 cannot do n>1

n = 5;
  
ellPIvalues22 = [-0.276 + 0.5554*i, -0.2714 + 0.5698*i, -0.2659 + 0.5854*i, -0.2591 + 0.6024*i, -0.2508 + 0.6209*i, -0.2404 + 0.6413*i, -0.2273 + 0.6638*i, -0.2105 + 0.6888*i, -0.1885 + 0.717*i, -0.1588 + 0.7488*i, -0.1165 + 0.7854*i, -0.0507 + 0.8279*i, -0.0799 + 0.0927*i, 0.2238 - 3.1*i, 0.1689 + 0.5812*i, 0.1417 + 0.5458*i, 0.1245 + 0.5222*i, 0.1124 + 0.5043*i, 0.1033 + 0.4899*i, 0.0961 + 0.4777*i, -0.0902 - 0.3456*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Incomplete PI(m<=1,b>pi/2,n>1): ')

if all(tryPI(m<=1)==ellPIvalues22(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues12(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to n>1\n')
end

fprintf('Incomplete PI(m<M, b>pi/2,n>1): ')

if all(tryPI(m<M)==ellPIvalues22(m<M)) 
  fprintf('passed\n')
  
elseif all(tryPI(m<M)==ellPIvalues22(m<M)) && all(tryPI(m==0)~=ellPIvalues22(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryPI(m<M))==real(ellPIvalues12(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed due to n>1\n')
end

fprintf('Incomplete PI(m>M, b>pi/2,n>1): ')

if all(tryPI(m>M)==ellPIvalues22(m>M)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>M))==real(ellPIvalues22(m>M)))
  fprintf('real components correct, missing complex\n')
else any(m>M)&& n>1;
  fprintf('failed due to a complex b input into elliptic3ic aswell as n>1\n')
end

%% incomplete Pi negative n wrong for all m

n = -5;

ellPIvalues23 = [0.5033, 0.513, 0.5236, 0.5355, 0.5489, 0.5642, 0.5821, 0.6036, 0.6304, 0.6656, 0.7167, 0.8085, inf, 0.7252 - 0.4*i, 0.6655 - 0.3738*i, 0.6275 - 0.3675*i, 0.5984 - 0.3655*i, 0.5744 - 0.3648*i, 0.5539 - 0.3644*i, 0.5359 - 0.3641*i, 0.5198 - 0.3637*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

fprintf('Incomplete PI(m<=1,b>pi/2,n<0): ')

if all(tryPI(m<=1)==ellPIvalues23(m<=1)) 
  fprintf('passed\n')
elseif all(real(tryPI(m<=1))==real(ellPIvalues23(m<=1)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m<M, b>pi/2,n<0): ')

if all(tryPI(m<M)==ellPIvalues23(m<M)) 
  fprintf('passed\n')
elseif all(tryPI(m<M)==ellPIvalues23(m<M)) && all(tryPI(m==0)~=ellPIvalues23(m==0))
  fprintf('problem at m=0 otherwise passed\n')
elseif all(real(tryPI(m<M))==real(ellPIvalues23(m<M)))
  fprintf('real components correct, missing complex\n')
else
  fprintf('failed\n')
end

fprintf('Incomplete PI(m>M, b>pi/2,n<0): ')

if all(tryPI(m>M)==ellPIvalues23(m>M)) 
  fprintf('passed\n')
elseif all(real(tryPI(m>M))==real(ellPIvalues23(m>M)))
  fprintf('real components correct, missing complex\n')
else any(m>M);
  fprintf('failed due to a complex b input into elliptic3ic\n')
end

