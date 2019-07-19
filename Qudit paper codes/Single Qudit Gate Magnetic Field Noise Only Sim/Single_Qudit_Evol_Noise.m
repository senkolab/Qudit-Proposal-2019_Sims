g0_factor=0.5; %Magnitude of hyperfine g-factor in Ba-137+ 6S level
t_ends=thetas./(2*pi*Rabi_frequency); %Gate time for each pulse in microseconds
t=0:deltat:sum(t_ends); %Initialize time array
t_noise=0:0.5*deltat:sum(t_ends); %Second time array with double the 
%resolution of t. To be used for Runge-Kutta-4 numerical integration.
psi00=zeros(8,1); %Initialize qudit state, including unencoded states.

uB=9.274009994*(10^(-24)); %Bohr magneton
Planck_h=6.62607015*(10^(-34)); %Planck constant

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

%Compute ideal Hamiltonian in time for all the pulses.
H_ideal=nan(numel(psi00),numel(psi00),1);
for pp=1:numel(thetas)
    [H_I,Off_resonant_freq]=Givens_Hamiltonian_TV(Rabi_frequency,spin_phases(pp),transitions(pp),g_array,m_array,off_resonant);
    if pp==1
        H_I=repmat(H_I,[1 1 numel(t_noise(t_noise>=0 & t_noise<=t_ends(pp)))]);
    else
        H_I=repmat(H_I,[1 1 numel(t_noise(t_noise>sum(t_ends(1:(pp-1))) & t_noise<=sum(t_ends(1:(pp)))))]);
    end
    H_ideal=cat(3,H_ideal,H_I);
end
H_ideal(:,:,1)=[];

%Get ideal time evolution operator
U_ideal=SingleQuditReconstruct(thetas,spin_phases,transitions,numel(psi00),0);

percent=0.1; %Initialize percent to track progress later in the code
tic
error_array=nan(n_iter,1); %Initialize the array of error figures 
%corresponding to the sample size
self_error_array=nan(n_iter,1); %Initialize the array of self-error, which 
%is the deviation from 1 when the magnitude of the output state is
%computed. Used to judge if time resolution is sufficient.

%Compute error figures
for hh=1:n_iter
    psi00=zeros(8,1);
    psi00(1)=1;
    psi_ideal=U_ideal*psi00;
    psi=nan(numel(psi00),numel(t));
    psi(:,1)=psi00;
    B_noise=Ornstein_Uhlenbeck(2000,2.7*(10^(-12)),t_noise); %Generate 
    %magnetic field noise
    B_noise=reshape(B_noise,[1 1 numel(B_noise)]);
    H_noise=2*pi.*repmat((10^(-6)).*(uB/Planck_h).*diag(g_array.*m_array),[1 1 numel(B_noise)])...
        .*repmat(B_noise,[numel(psi00) numel(psi00) 1]); %Compute noise Hamiltonian
    H_net=H_ideal+H_noise; %Net Hamiltonian
for h=2:numel(t) %Runge-Kutta-4 numerical integration
    k1=-1i.*(H_net(:,:,2*(h-1)-1)*psi(:,h-1)).*(t(h)-t(h-1));
    k2=-1i.*(H_net(:,:,2*(h-1))*(psi(:,h-1)+k1./2)).*(t(h)-t(h-1));
    k3=-1i.*(H_net(:,:,2*(h-1))*(psi(:,h-1)+k2./2)).*(t(h)-t(h-1));
    k4=-1i.*(H_net(:,:,2*(h-1)+1)*(psi(:,h-1)+k3)).*(t(h)-t(h-1));
    psi(:,h)=psi(:,h-1)+(1/6).*(k1+2.*k2+2.*k3+k4);
end
error_array(hh)=1-abs(psi_ideal'*psi(:,end))^2;
self_error_array(hh)=1-psi(:,end)'*psi(:,end);
if hh==round(percent*n_iter)
    display([num2str(100*percent) '% done']);
    percent=percent+0.1;
end
end
toc