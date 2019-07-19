function [H_I,Off_resonant_freq]=Givens_Hamiltonian_TV(Rabi_frequency,spin_phase,transition,g_array,m_array,off_resonant)
%Input definition:
%Rabi_frequency: Rabi frequency of desired transition in MHz. Scalar
%spin_phase: Spin phase of the desired transition in radian. Scalar
%transition: Selects the transition from (transition-1) to (transition).
%Integer
%g_array: g-factors corresponding to the quantum state. 1D array.
%m_array: m numbers corresponding to the quantum state. 1D array.
%off_resonant: sets output to give Hamiltonian without off-resonance if set
%to 0. Code outputs Hamiltonian with off-resonance otherwise.

%Output definition:
%H_I: Rabi frequencies of the Hamiltonian in 2*pi MHz. 2D matrix
%Off_resonant_freq: Off-resonant frequencies in 2*pi MHz. 2D matrix

%Initialize output matrices
H_I=zeros(numel(g_array),numel(g_array));
Rabi_Coeff=nan(numel(g_array),numel(g_array));

%Get Clebsch-Gordan matrix
[A,J_index,Jm_index,M1_index,M2_index]=GetClebschGordan(1,1);

%Compute Clebsch-Gordan coefficients for off-resonant transitions
for h=1:size(Rabi_Coeff,1)
    for hh=1:size(Rabi_Coeff,2)
        if abs(m_array(h)-m_array(hh))<=1 && g_array(h)>g_array(hh)
            C1=A(M1_index==m_array(hh) & M2_index==(m_array(h)-m_array(hh)),J_index==2 & Jm_index==m_array(h));
        elseif abs(m_array(h)-m_array(hh))<=1 && g_array(h)<g_array(hh)
            C1=A(M1_index==m_array(h) & M2_index==(m_array(hh)-m_array(h)),J_index==2 & Jm_index==m_array(hh));
        else
            C1=0;
        end
        if isempty(C1)
            C1=0;
        end
        Rabi_Coeff(h,hh)=C1;
    end
end

%Compute H_I
%Select transitions allowed by transition rule
for h=1:numel(g_array)
    H_dummy=zeros(1,numel(g_array));
    H_dummy(g_array==-g_array(h) & abs(m_array-m_array(h))<=1)=Rabi_frequency/2;
    H_I(h,:)=H_dummy;
end
%Incorporate spin phase
H_I=triu(H_I).*exp(-1i*spin_phase)+tril(H_I).*exp(1i*spin_phase);
%Compute Rabi frequencies for off-resonant transitions
Rabi_Coeff=Rabi_Coeff./Rabi_Coeff(transition,transition+1);
H_I=Rabi_Coeff.*H_I;
H_I=2.*pi.*H_I;

%Compute off-resonant frequencies
[g_sign_array_column,g_sign_array_row]=meshgrid(sign(g_array),sign(g_array));
[m_array_column,m_array_row]=meshgrid(m_array,m_array);
Off_resonant_freq=g_sign_array_row.*m_array_row-g_sign_array_column.*m_array_column;
Multiple_ref=g_sign_array_row(transition)*Off_resonant_freq(transition,transition+1);
Off_resonant_freq=(Off_resonant_freq-g_sign_array_row.*Multiple_ref).*3.2875;
Off_resonant_freq(g_sign_array_row==g_sign_array_column)=0;
Off_resonant_freq=2.*pi.*Off_resonant_freq;

%Compute ideal Hamiltonian if off_resonant is set to 0
if off_resonant==0
    H_I_dummy=H_I(transition,transition+1);
    H_I_dummy_2=H_I(transition+1,transition);
    H_I=zeros(numel(g_array),numel(g_array));
    H_I(transition,transition+1)=H_I_dummy;
    H_I(transition+1,transition)=H_I_dummy_2;
    Off_resonant_freq=zeros(size(H_I));
end