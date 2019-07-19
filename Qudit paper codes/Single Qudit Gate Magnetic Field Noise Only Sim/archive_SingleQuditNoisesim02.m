qudit_level=3;
g0_factor=0.5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/4 pi/2 pi/2 atan(sqrt(2)) pi/4]; %Fourier Transform 3-level
spin_phases=-[pi/2 3*pi/2 5*pi/6 pi/2 2*pi/3 7*pi/6 7*pi/6]; %Fourier Transform 3-level
transitions=[1 1 1 2 2 2 1]; %Fourier Transform 3-level
deltat=0.0002;
n_iter=300;
off_resonant=0;

Single_Qudit_Evol_Noise;

error_array_F3_RealNoise_on_Rabi_0o01=error_array;
self_error_array_F3_RealNoise_on_Rabi_0o01=self_error_array;

qudit_level=5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/4 pi/2 pi/2 0.95532 0.60641 pi/2 pi/2 pi/3 0.85289 0.60641 pi/2 pi/2 1.10714 pi/3 0.95532 pi/4]; %Fourier Transform 5-level
spin_phases=-[pi/2 3.30265 0.63627 pi/2 6.18626 1.53005 4.57966 pi/2 pi/2 1.981884 3.74954 3.69336 pi/2 9*pi/10 9*pi/10 9*pi/10 9*pi/10 9*pi/10]; %Fourier Transform 5-level
transitions=[1 1 1 2 2 2 1 3 3 3 2 1 4 4 4 3 2 1]; %Fourier Transform 5-level
n_iter=300;

Single_Qudit_Evol_Noise;

error_array_F5_RealNoise_on_Rabi_0o01=error_array;
self_error_array_F5_RealNoise_on_Rabi_0o01=self_error_array;