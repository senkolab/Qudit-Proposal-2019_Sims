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
n_average_COM=0.1; %Average COM phonon number
n_cycle=1; %Number of cycles in phase space, typically denoted K in literatures
heattime=0; %The fraction of gate time when a phonon jump due to motional heating occurs
mag_shift_on=1; 

%Fidelity with all error sources except phonon heating
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_real=fidel;

%Fidelity without magnetic field shift
mag_shift_on=0; 
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_no_mag_shift=fidel;
mag_shift_on=1;

%Fidelity with large detuning to minimize error from RWA
fq_COM=50;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_noRWA1=fidel;

%Fidelity with larger detuning to ensure error from RWA is eliminated to
%4th significant figure
fq_COM=60;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_noRWA2=fidel;

%Fidelity with pure input state
fq_COM=2;
fq_Tilt=fq_COM-0.2;
n_average_COM=0;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_nothermal=fidel;

%Fidelity without tilt mode
n_average_COM=0.1;
MS_gate_qudit_2modes_ode113_COM_lowmem;
fidel_5qudit_nospec=fidel;

%Fidelity without error from Lamb-Dicke approximation
MS_gate_qudit_2modes_ode113_LDA_lowmem;
fidel_5qudit_noLD=fidel;


mass=0.137./(6.022140857*(10^23)); %Mass for Ba-137
lambda_Raman=532; %Raman laser wavelength
hbar=1.0545718*10^(-34); %Plank's constant

qudit_level=5; %qudit level
gate_time=100; %Gate time in micro-seconds
fq_COM=2; %COM motional frequency in MHz
fq_Tilt=fq_COM-0.2; %Tilt mode motional frequency in MHz
eta_COM=sqrt(2)*(2*pi/(lambda_Raman*(10^(-9))))*sqrt(hbar/(2*2*pi*fq_COM*(10^6)*2*mass)); %COM Lamb-Dicke parameter
n_average_COM=0.1; %Average COM phonon number
n_cycle=1; %Number of cycles in phase space, typically denoted K in literatures
heattime=0.5; %The fraction of gate time when a phonon jump due to motional heating occurs
mag_shift_on=1; 

%Fidelity with all error sources except phonon heating
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_real_heat=fidel;

%Fidelity without magnetic field shift
mag_shift_on=0; 
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_no_mag_shift_heat=fidel;
mag_shift_on=1;

%Fidelity with large detuning to minimize error from RWA
fq_COM=50;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_noRWA1_heat=fidel;

%Fidelity with larger detuning to ensure error from RWA is eliminated to
%4th significant figure
fq_COM=60;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_noRWA2_heat=fidel;

%Fidelity with pure input state
fq_COM=2;
fq_Tilt=fq_COM-0.2;
n_average_COM=0;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_5qudit_nothermal_heat=fidel;

%Fidelity without tilt mode
n_average_COM=0.1;
MS_gate_qudit_2modes_ode113_COM_lowmem;
fidel_5qudit_nospec_heat=fidel;

%Fidelity without error from Lamb-Dicke approximation
MS_gate_qudit_2modes_ode113_LDA_lowmem;
fidel_5qudit_noLD_heat=fidel;

disp('5-level done')

mass=0.137./(6.022140857*(10^23)); %Mass for Ba-137
lambda_Raman=532; %Raman laser wavelength
hbar=1.0545718*10^(-34); %Plank's constant

qudit_level=3; %qudit level
gate_time=100; %Gate time in micro-seconds
fq_COM=2; %COM motional frequency in MHz
fq_Tilt=fq_COM-0.2; %Tilt mode motional frequency in MHz
eta_COM=sqrt(2)*(2*pi/(lambda_Raman*(10^(-9))))*sqrt(hbar/(2*2*pi*fq_COM*(10^6)*2*mass)); %COM Lamb-Dicke parameter
n_average_COM=0.1; %Average COM phonon number
n_cycle=1; %Number of cycles in phase space, typically denoted K in literatures
heattime=0; %The fraction of gate time when a phonon jump due to motional heating occurs
mag_shift_on=1; 

%Fidelity with all error sources except phonon heating
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_real=fidel;

%Fidelity without magnetic field shift
mag_shift_on=0; 
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_no_mag_shift=fidel;
mag_shift_on=1;

%Fidelity with large detuning to minimize error from RWA
fq_COM=50;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_noRWA1=fidel;

%Fidelity with larger detuning to ensure error from RWA is eliminated to
%4th significant figure
fq_COM=60;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_noRWA2=fidel;

%Fidelity with pure input state
fq_COM=2;
fq_Tilt=fq_COM-0.2;
n_average_COM=0;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_nothermal=fidel;

%Fidelity without tilt mode
n_average_COM=0.1;
MS_gate_qudit_2modes_ode113_COM_lowmem;
fidel_3qudit_nospec=fidel;

%Fidelity without error from Lamb-Dicke approximation
MS_gate_qudit_2modes_ode113_LDA_lowmem;
fidel_3qudit_noLD=fidel;


mass=0.137./(6.022140857*(10^23)); %Mass for Ba-137
lambda_Raman=532; %Raman laser wavelength
hbar=1.0545718*10^(-34); %Plank's constant

qudit_level=3; %qudit level
gate_time=100; %Gate time in micro-seconds
fq_COM=2; %COM motional frequency in MHz
fq_Tilt=fq_COM-0.2; %Tilt mode motional frequency in MHz
eta_COM=sqrt(2)*(2*pi/(lambda_Raman*(10^(-9))))*sqrt(hbar/(2*2*pi*fq_COM*(10^6)*2*mass)); %COM Lamb-Dicke parameter
n_average_COM=0.1; %Average COM phonon number
n_cycle=1; %Number of cycles in phase space, typically denoted K in literatures
heattime=0.5; %The fraction of gate time when a phonon jump due to motional heating occurs
mag_shift_on=1; 

%Fidelity with all error sources except phonon heating
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_real_heat=fidel;

%Fidelity without magnetic field shift
mag_shift_on=0; 
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_no_mag_shift_heat=fidel;
mag_shift_on=1;

%Fidelity with large detuning to minimize error from RWA
fq_COM=50;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_noRWA1_heat=fidel;

%Fidelity with larger detuning to ensure error from RWA is eliminated to
%4th significant figure
fq_COM=60;
fq_Tilt=fq_COM-0.2;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_noRWA2_heat=fidel;

%Fidelity with pure input state
fq_COM=2;
fq_Tilt=fq_COM-0.2;
n_average_COM=0;
MS_gate_qudit_2modes_ode113_mag_shift_lowmem;
fidel_3qudit_nothermal_heat=fidel;

%Fidelity without tilt mode
n_average_COM=0.1;
MS_gate_qudit_2modes_ode113_COM_lowmem;
fidel_3qudit_nospec_heat=fidel;

%Fidelity without error from Lamb-Dicke approximation
MS_gate_qudit_2modes_ode113_LDA_lowmem;
fidel_3qudit_noLD_heat=fidel;