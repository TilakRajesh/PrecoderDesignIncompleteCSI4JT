function N0 = CoMP_Rx_NoiseInit()
T = 290; %Operating temperature , in degrees Kelvin, T = N/kB
kB = 1.38e-23;
B = 10e6; % noise bandwidth
N0 = T*kB*B; % lineal 
end