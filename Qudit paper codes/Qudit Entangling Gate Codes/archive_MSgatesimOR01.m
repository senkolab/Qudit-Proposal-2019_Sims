options=odeset('RelTol',1e-5);

%Repeat for 5-level qudit
mass=0.137./(6.022140857*(10^23)); %Mass for Ba-137
lambda_Raman=532; %Raman laser wavelength
hbar=1.0545718*10^(-34); %Plank's constant

qudit_level=5; %qudit level
gate_time=100; %Gate time in micro-seconds
fq_COM=2; %COM motional frequency in MHz
fq_Tilt=fq_COM-0.2; %Tilt mode motional frequency in MHz
eta_COM=sqrt(2)*(2*pi/(lambda_Raman*(10^(-9))))*sqrt(hbar/(2*2*pi*fq_COM*(10^6)*2*mass)); %COM Lamb-Dicke parameter
n_average_COM=0; %Average COM phonon number
n_cycle=1; %Number of cycles in phase space, typically denoted K in literatures
heattime=0; %The fraction of gate time when a phonon jump due to motional heating occurs
mag_shift_on=1; 
w_split=3.2875*2*pi;
Cfactor=[3 1/3 1/3 3];

%Fidelity with all error sources except phonon heating
MS_gate_qudit_2modes_ode113_ideal_OR_lowmem_2;
fidel_5qudit_OR=fidel;