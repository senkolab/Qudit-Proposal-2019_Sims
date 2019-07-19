uB=9.274009994*(10^(-24)); %Bohr magneton
Planck_h=6.62607015*(10^(-34)); %Planck constant

t_ends=thetas./(2*pi*Rabi_frequency); %Gate time for each pulse
deltat=0.00001; %Unused
error_array=nan(n_iter,1); %Initialize array of error values of the sample 
%size

%Initialize g-factor and m numbers corresponding to encoding schemes for 3
%or 5-level qudits. The first d elements in m_array and g_array correspond
%to the encoded states. The rest are the unencoded states.
if qudit_level==5
    m_array=[-2 -1 0 1 2 -1 0 1];
    g_array=g0_factor.*[1 -1 1 -1 1 1 -1 1];
elseif qudit_level==3
    m_array=[-1 0 1 -2 -1 0 1 2];
    g_array=[1 -1 1 1 -1 1 -1 1];
end
H_noise_constant=2*pi.*(10^(-6)).*(uB/Planck_h).*diag(g_array.*m_array)...
    .*2.7.*(10^(-12)); %Generate the Hamiltonian corresponding to a 
%magnetic field offset
percent=0.1; %Initialize percent = 10% to track progress of simulations
display('Initialized')
tic

for p=1:n_iter
psi00=zeros(8,1); %Initialize input state
psi00(1:qudit_level)=[rand(qudit_level,1)+1i.*rand(qudit_level,1)]; 
%Randomize input state
psi00=psi00./sqrt(psi00'*psi00); %Normalize input state
psi0=psi00; %Redundant reasignment
t0=0; %Start at t = 0 for first pulse
for pp=1:numel(thetas)
    [H_I,Off_resonant_freq]=Givens_Hamiltonian_TV(Rabi_frequency,spin_phases(pp),transitions(pp),g_array,m_array,off_resonant);
    %Get Rabi frequencies of Hamiltonian in H_I and the off-resonant
    %frequencies in Off_resonant_freq.
    [tm,psi_t]=ode113(@(t,psi) Qudit_Givens(t,psi,H_I+H_noise_on.*H_noise_constant,Off_resonant_freq),...
        [t0 t0+t_ends(pp)],psi0,options); %Compute time evolution of state
    psi_t=psi_t.';
    psi0=psi_t(:,end); %Set input state for next pulse
    t0=tm(end); %Set start time for next pulse
end

U_ideal=SingleQuditReconstruct(thetas,spin_phases,transitions,numel(psi00),0);
%Get ideal time evolution operator
psi_ideal=U_ideal*psi00; %Get ideal output state
error_array(p)=1-abs(psi_ideal'*psi_t(:,end))^2; %Compute error
if p==round(percent*n_iter) %Track progress
    display([num2str(100*percent) '% done']);
    percent=percent+0.1;
    toc
end
end