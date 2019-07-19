options=odeset('AbsTol',1e-10,'RelTol',1e-10);
qudit_level=3;
g0_factor=0.5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/4 pi/2 pi/2 atan(sqrt(2)) pi/4]; %Fourier Transform 3-level
spin_phases=-[pi/2 3*pi/2 5*pi/6 pi/2 2*pi/3 7*pi/6 7*pi/6]; %Fourier Transform 3-level
transitions=[1 1 1 2 2 2 1]; %Fourier Transform 3-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_F3_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_F3_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_F3_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_F3_Noise_off_Rabi_0o1=error_array;

qudit_level=5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/4 pi/2 pi/2 0.95532 0.60641 pi/2 pi/2 pi/3 0.85289 0.60641 pi/2 pi/2 1.10714 pi/3 0.95532 pi/4]; %Fourier Transform 5-level
spin_phases=-[pi/2 3.30265 0.63627 pi/2 6.18626 1.53005 4.57966 pi/2 pi/2 1.981884 3.74954 3.69336 pi/2 9*pi/10 9*pi/10 9*pi/10 9*pi/10 9*pi/10]; %Fourier Transform 5-level
transitions=[1 1 1 2 2 2 1 3 3 3 2 1 4 4 4 3 2 1]; %Fourier Transform 5-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_F5_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_F5_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_F5_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_F5_Noise_off_Rabi_0o1=error_array;

qudit_level=3;
g0_factor=0.5;
Rabi_frequency=0.01;
thetas=2.*[pi pi/2 pi/2]; %X 3-level
spin_phases=-[0 pi/2 pi/2]; %X 3-level
transitions=[1 2 1]; %X 3-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_X3_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_X3_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_X3_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_X3_Noise_off_Rabi_0o1=error_array;

qudit_level=3;
g0_factor=0.5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/2 pi/2]; %Y 3-level
spin_phases=-[pi/2 7*pi/6 pi/2 pi/2]; %Y 3-level
transitions=[1 1 2 1]; %Y 3-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_Y3_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Y3_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_Y3_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Y3_Noise_off_Rabi_0o1=error_array;

qudit_level=3;
g0_factor=0.5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2]; %Z 3-level
spin_phases=-[pi/2 pi/6]; %Z 3-level
transitions=[2 2]; %Z 3-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_Z3_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Z3_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_Z3_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Z3_Noise_off_Rabi_0o1=error_array;

qudit_level=3;
g0_factor=0.5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2]; %T 3-level
spin_phases=-[pi/2 31*pi/18]; %T 3-level
transitions=[2 2]; %T 3-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_T3_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_T3_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_T3_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_T3_Noise_off_Rabi_0o1=error_array;

qudit_level=5;
Rabi_frequency=0.01;
thetas=2.*[pi pi pi/2 pi/2 pi/2 pi/2]; %X 5-level
spin_phases=-[0 0 pi/2 pi/2 pi/2 pi/2]; %X 5-level
transitions=[1 3 4 3 2 1]; %X 5-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_X5_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_X5_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_X5_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_X5_Noise_off_Rabi_0o1=error_array;

qudit_level=5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %Y 5-level
spin_phases=-[pi/2 9*pi/10 pi/2 7*pi/10 pi/2 9*pi/10 pi/2 pi/2 pi/2 pi/2]; %Y 5-level
transitions=[1 1 2 2 3 3 4 3 2 1]; %Y 5-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_Y5_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Y5_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_Y5_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Y5_Noise_off_Rabi_0o1=error_array;

qudit_level=5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %Z 5-level
spin_phases=-[pi/2 19*pi/10 pi/2 7*pi/10 pi/2 19*pi/10]; %Z 5-level
transitions=[2 2 3 3 4 4]; %Z 5-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_Z5_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Z5_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_Z5_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_Z5_Noise_off_Rabi_0o1=error_array;

qudit_level=5;
Rabi_frequency=0.01;
thetas=2.*[pi/2 pi/2 pi/2 pi/2 pi/2 pi/2]; %T 5-level
spin_phases=-[pi/2 7*pi/10 pi/2 3*pi/10 pi/2 11*pi/10]; %T 5-level
transitions=[2 2 3 3 4 4]; %T 5-level
H_noise_on=1;
n_iter=500;
off_resonant=1;

Single_Qudit_Off_Sim;
error_array_T5_Noise_on_Rabi_0o01=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_T5_Noise_off_Rabi_0o01=error_array;

H_noise_on=1;
Rabi_frequency=0.1;
Single_Qudit_Off_Sim;
error_array_T5_Noise_on_Rabi_0o1=error_array;

H_noise_on=0;
Single_Qudit_Off_Sim;
error_array_T5_Noise_off_Rabi_0o1=error_array;