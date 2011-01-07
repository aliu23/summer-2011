%% Complete F (true for all) 

m = -5:0.5:5;
ellKvalues = [0.9555, 0.981, 1.0095, 1.0416, 1.0783, 1.1209, 1.1714, 1.233, 1.311, 1.4157, 1.5708, 1.8541, inf, 1.6566 - 1.4157*i, 1.311 - 1.311*i, 1.1242 - 1.233*i, 1.0011 - 1.1714*i, 0.9117 - 1.1209*i, 0.8429 - 1.0783*i, 0.7877 - 1.0416*i, 0.7422 - 1.0095*i];

tryk = nan(size(m));

for ii = 1:length(m)

    tryk(ii) = elliptic12(m(ii));
    
end

tryk = round(10000*tryk)/10000;

tryk==ellKvalues % creates a vector of true's and falses

%% Complete E (true for all)

m = -5:0.5:5;
ellEvalues = [2.8302, 2.7347, 2.6352, 2.5312, 2.4221, 2.3069, 2.1844, 2.053, 1.9101, 1.7518, 1.5708, 1.3506, 1., 0.7163 + 0.336*i, 0.5991 + 0.5991*i, 0.5263 + 0.82*i, 0.4752 + 1.013*i, 0.4367 + 1.186*i, 0.4063 + 1.3439*i, 0.3815 + 1.4897*i, 0.3608 + 1.6257*i];

tryE = nan(size(m));

for ii = 1:length(m)

    [~, tryE(ii)] = elliptic12(m(ii));
    
end

tryE = round(10000*tryE)/10000;

tryE==ellEvalues

%% complete pi ,n<1, input must be real, elliptic3 cannot have complex inputs
%we need elliptic3i

n = 0.5;
m = -5:0.5:5;

ellPIvalues1 = [1.2566, 1.2943, 1.3367, 1.3847, 1.44, 1.5048, 1.5823, 1.6779, 1.8005, 1.9679, 2.2214, 2.7013, inf, 2.0855 - 2.4344*i, 1.5327 - 2.1041*i, 1.2659 - 1.8989*i, 1.1017 - 1.754*i, 0.9879 - 1.6439*i, 0.9032 - 1.5561*i, 0.837 - 1.4838*i, 0.7835 - 1.4227*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
        tryPI(ii) = elliptic3(m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues1
%% Complete Pi n>1 2 problems no elliptic3i and can only output real values
%for real inputs

n = 5
m = -5:0.5:5;

ellPIvalues2 = [0.2271 - 0.5554*i, 0.2204 - 0.5698*i, 0.2124 - 0.5854*i, 0.2028 - 0.6024*i, 0.1911 - 0.6209*i, 0.1766 - 0.6413*i, 0.1585 - 0.6638*i, 0.1354 - 0.6888*i, 0.1048 - 0.717*i, 0.0626 - 0.7488*i, -0.7854*i, -0.1092 - 0.8279*i, inf, -0.2238 - 0.469*i, -0.1689 - 0.4558*i, -0.1417 - 0.4448*i, -0.1245 - 0.4353*i, -0.1124 - 0.4269*i, -0.1033 - 0.4194*i, -0.0961 - 0.4126*i, -0.0902 - 0.4064*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
        tryPI(ii) = elliptic3(m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI==ellPIvalues2;

real(tryPI)==real(ellPIvalues2);

tryPI

%% complete Pi negative n (no elliptic3i prob)

n = -5;
m = -5:0.5:5;

ellPIvalues3 = [0.4717, 0.48, 0.489, 0.499, 0.5102, 0.523, 0.5377, 0.555, 0.5762, 0.6033, 0.6413, 0.7049, inf, 0.7252 - 0.2799*i, 0.6655 - 0.2928*i, 0.6275 - 0.3023*i, 0.5984 - 0.3094*i, 0.5744 - 0.3148*i, 0.5539 - 0.3189*i, 0.5359 - 0.322*i, 0.5198 - 0.3244*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues3
%% incomplete F elliptic12i cannot output the imaginary part for m values
% 4 and over,, real part ok

 b = 0.5
 m = -5:0.5:5
ellFvalues1 = [0.4325, 0.4374, 0.4426, 0.4481, 0.454, 0.4602, 0.4669, 0.4742, 0.482, 0.4906, 0.5, 0.5105, 0.5222, 0.5357, 0.5514, 0.5702, 0.5938, 0.6258, 0.6774, 0.7877 - 0.0986*i, 0.7422 - 0.1898*i];

tryF = nan(size(m));

for ii = 1:length(m)

    tryF(ii) = elliptic12(b,m(ii));
    
end

tryF = round(10000*tryF)/10000;

tryF==ellFvalues1

real(tryF)==real(ellFvalues1);

imag(tryF)==imag(ellFvalues1);
%% incomplete E (same prob as F)

b=0.5;
m = -5:0.5:5;
ellEvalues1 = [0.5864, 0.5787, 0.5707, 0.5627, 0.5544, 0.5459, 0.5372, 0.5283, 0.5192, 0.5097, 0.5, 0.4899, 0.4794, 0.4685, 0.457, 0.4448, 0.4319, 0.4177, 0.4018, 0.3815 + 0.0011*i, 0.3608 + 0.0093*i];

tryE = nan(size(m));

for ii = 1:length(m)

    [~, tryE(ii)] = elliptic12(b,m(ii));
    
end

tryE = round(10000*tryE)/10000;

tryE

tryE==ellEvalues1

%% incomplete Pi n<1 (elliptic 3 cannot imput complex no.)

b=0.5
n = 0.5
m = -5:0.5:5

ellPIvalues11 = [0.449, 0.4543, 0.4598, 0.4657, 0.4719, 0.4786, 0.4858, 0.4936, 0.502, 0.5112, 0.5213, 0.5326, 0.5453, 0.5598, 0.5768, 0.5973, 0.623, 0.658, 0.7149, 0.837 - 0.1111*i, 0.7835 - 0.2121*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues11


%% Pi when n>1 Cant do any as n>1
b=0.5
n = 5
m = -5:0.5:5

ellPIvalues12 = [0.6705 - 0.5554*i, 0.679 - 0.5698*i, 0.688 - 0.5854*i, 0.6974 - 0.6024*i, 0.7073 - 0.6209*i, 0.7178 - 0.6413*i, 0.7288 - 0.6638*i, 0.7405 - 0.6888*i, 0.7528 - 0.717*i, 0.7658 - 0.7488*i, 0.7795 - 0.7854*i, 0.7937 - 0.8279*i, 0.8084 + 0.8781*i, 0.823 + 2.8162*i, 0.8365 + 3.0418*i, 0.8466 + 3.3322*i, 0.847 + 3.7255*i, 0.8188 + 4.3018*i, 0.6791 + 5.2686*i, -0.0961 - 1.6804*i, -0.0902 - 1.3637*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues12

%% incomplete F b=2 (for all m<0 the output is wrong ?) 
%for m>1 correct up until we encounter complex outputs then it just spits
%out the reals


b=2
m = -5:0.5:5;
ellFvalues2 = [1.1354, 1.1688, 1.2063, 1.2489, 1.2979, 1.3554, 1.4244, 1.5095, 1.6192, 1.7697, 2., 2.4444, inf, 1.6566 - 2.0956*i, 1.311 - 1.7707*i, 1.1242 - 1.6035*i, 1.0011 - 1.4903*i, 0.9117 - 1.405*i, 0.8429 - 1.337*i, 0.7877 - 1.2807*i, 0.7422 - 1.2329*i];
tryF = nan(size(m));

for ii = 1:length(m)

    tryF(ii) = elliptic12(b,m(ii));
    
end

tryF = round(10000*tryF)/10000;

tryF

tryF==ellFvalues2

real(tryF)==real(ellFvalues2);

imag(tryF)==imag(ellFvalues2);

%% incomplete E doesnt work for any m<0??

b=2;
m = -5:0.5:5;
ellEvalues2 = [3.855, 3.7163, 3.5717, 3.4203, 3.2611, 3.0926, 2.9129, 2.7194, 2.508, 2.2722, 2., 1.6629, 1.0907, 0.7163 + 0.6099*i, 0.5991 + 1.0013*i, 0.5263 + 1.3184*i, 0.4752 + 1.5919*i, 0.4367 + 1.8354*i, 0.4063 + 2.0568*i, 0.3815 + 2.261*i, 0.3608 + 2.4513*i];

tryE = nan(size(m));

for ii = 1:length(m)

    [~, tryE(ii)] = elliptic12(b,m(ii));
    
end

tryE = round(10000*tryE)/10000;

tryE

tryE==ellEvalues2

real(tryE)==real(ellEvalues2);

imag(tryE)==imag(ellEvalues2);

%% incomplete Pi n<1 Doesnt work for all m same no elliptic 3i problem
b=2;
n = 0.5;
m = -5:0.5:5;
ellPIvalues21 = [1.5967, 1.6493, 1.7087, 1.7767, 1.8554, 1.9483, 2.0607, 2.2008, 2.3834, 2.6376, 3.0338, 3.8199, inf, 2.0855 - 3.7147*i, 1.5327 - 2.9716*i, 1.2659 - 2.5984*i, 1.1017 - 2.3562*i, 0.9879 - 2.1806*i, 0.9032 - 2.045*i, 0.837 - 1.9357*i, 0.7835 - 1.8449*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues21

%% incomplete Pi n>1 cannot do n>1

  n = 5;
  m = -5:0.5:5;
ellPIvalues22 = [-0.276 + 0.5554*i, -0.2714 + 0.5698*i, -0.2659 + 0.5854*i, -0.2591 + 0.6024*i, -0.2508 + 0.6209*i, -0.2404 + 0.6413*i, -0.2273 + 0.6638*i, -0.2105 + 0.6888*i, -0.1885 + 0.717*i, -0.1588 + 0.7488*i, -0.1165 + 0.7854*i, -0.0507 + 0.8279*i, -0.0799 + 0.0927*i, 0.2238 - 3.1*i, 0.1689 + 0.5812*i, 0.1417 + 0.5458*i, 0.1245 + 0.5222*i, 0.1124 + 0.5043*i, 0.1033 + 0.4899*i, 0.0961 + 0.4777*i, -0.0902 - 0.3456*i];

%% incomplete Pi negative n wrong for all m
b=2
n = -5
m = -5:0.5:5
ellPIvalues23 = [0.5033, 0.513, 0.5236, 0.5355, 0.5489, 0.5642, 0.5821, 0.6036, 0.6304, 0.6656, 0.7167, 0.8085, inf, 0.7252 - 0.4*i, 0.6655 - 0.3738*i, 0.6275 - 0.3675*i, 0.5984 - 0.3655*i, 0.5744 - 0.3648*i, 0.5539 - 0.3644*i, 0.5359 - 0.3641*i, 0.5198 - 0.3637*i];

tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues23
%% incomplete F (b not between 0 and pi/2)

b=5;
m = -5:0.5:5;
ellFvalues3 = [2.9853, 3.067, 3.1584, 3.2617, 3.3801, 3.5179, 3.6819, 3.8825, 4.1379, 4.4831, 5., 5.9636, inf, 4.9699 - 4.6726*i, 3.9331 - 4.2291*i, 3.3726 - 3.9395*i, 3.0032 - 3.722*i, 2.7351 - 3.5482*i, 2.5286 - 3.404*i, 2.3632 - 3.2812*i, 2.2266 - 3.1747*i];

tryF = nan(size(m));

for ii = 1:length(m)

    tryF(ii) = elliptic12(b,m(ii));
    
end

tryF = round(10000*tryF)/10000;

tryF

tryF==ellFvalues3;

real(tryF)==real(ellFvalues3)

imag(tryF)==imag(ellFvalues3);


%% incomplete E  even the real part is wrong for m>1 cases, m<0 is
%completely wrong (as phase should be between 0 and pi/2)
b=5;
m = -5:0.5:5;

ellEvalues3 = [9.1871, 8.8709, 8.5416, 8.1973, 7.8357, 7.4536, 7.0469, 6.61, 6.1343, 5.606, 5., 4.258, 3.0411, 2.1489 + 1.2029*i, 1.7972 + 2.0768*i, 1.579 + 2.8041*i, 1.4257 + 3.4374*i, 1.3101 + 4.0041*i, 1.2189 + 4.5206*i, 1.1445 + 4.9976*i, 1.0823 + 5.4425*i];

tryE = nan(size(m));

for ii = 1:length(m)

    [~, tryE(ii)] = elliptic12(b,m(ii));
    
end

tryE = round(10000*tryE)/10000;

tryE

tryE==ellEvalues3

real(tryE)==real(ellEvalues3);

imag(tryE)==imag(ellEvalues3);

%% incomplete pi (cannot evaluate b>pi/2)
b=5;
n = 0.5;
m = -5:0.5:5;
ellPIvalues31 = [4.0013, 4.1246, 4.2633, 4.4211, 4.6031, 4.8169, 5.0734, 5.3908, 5.8003, 6.3631, 7.2247, 8.8861, inf, 6.2566 - 8.1313*i, 4.5982 - 6.8888*i, 3.7976 - 6.165*i, 3.3051 - 5.6665*i, 2.9638 - 5.2929*i, 2.7097 - 4.9979*i, 2.5111 - 4.7563*i, 2.3505 - 4.5532*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues31

%% incomplete Pi
b=5;
n = 5;
m = -5:0.5:5;
ellPIvalues32 = [-0.2579 + 0.5554*i, -0.2525 + 0.5698*i, -0.2461 + 0.5854*i, -0.2383 + 0.6024*i, -0.2287 + 0.6209*i, -0.2169 + 0.6413*i, -0.202 + 0.6638*i, -0.1829 + 0.6888*i, -0.1578 + 0.717*i, -0.1237 + 0.7488*i, -0.0745 + 0.7854*i, 0.0052 + 0.8279*i, 0.2006 + 0.8781*i, 0.2238 - 3.1756*i, 0.1689 + 0.5326*i, 0.1417 + 0.5072*i, 0.1245 + 0.4892*i, 0.1124 + 0.475*i, 0.1033 + 0.4633*i, 0.0961 + 0.4532*i, -0.2706 - 1.1814*i];
tryPI = nan(size(m));

for ii = 1:length(m)
    try
    tryPI(ii) = elliptic3(b,m(ii),n);
    end
end

tryPI = round(10000*tryPI)/10000;

tryPI

tryPI==ellPIvalues32

% this is the end
