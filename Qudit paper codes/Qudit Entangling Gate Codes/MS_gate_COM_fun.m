function dpsidt = MS_gate_COM_fun(t,psi,qudit_level,gate_time,n_cycle,eta_COM,w_COM,Omega_bar,...
    State_raise,State_lower,State_I,Motion_raise_COM,Motion_lower_COM,...
    Motion_I_Tilt)

%Input definition
%t: time in micro-seconds, scalar
%psi: qudit state, 1D array
%gate_time: gate_time, scalar in micro-seconds
%n_cycle: number of cycles in phase space, usually denoted as K in
%literatures, scalar
%eta_COM: COM Lamb-Dicke parameter, scalar
%w_COM: COM motional angular frequency in MHz, scalar
%w_Tilt: Tilt motional angular frequency in MHz, scalar
%Omega_bar: Rabi frequency in MHz, scalar
%State_raise: Raising operator with appropriate Rabi frequency coefficients
%for the internal state, 2D matrix
%State_lower: Lowering operator with appropriate Rabi frequency coefficients
%for the internal state, 2D matrix
%State_I: Identity operator for the internal state, 2D matrix
%Motion_raise_COM: Raising operator for the COM Fock state, 2D matrix
%Motion_lower_COM: Lowering operator for the COM Fock state, 2D matrix
%Motion_I_Tilt: Identity operator for the tilt mode Fock state, 2D matrix

%Output definition
%dpsidt: Derivative of state psi, 1D array

delta=w_COM+n_cycle*2*pi/gate_time; %Laser detuning frequency in MHz
%eta_Tilt=eta_COM*sqrt(w_COM/w_Tilt);

%Define operators to obtain the Hamiltonian as outlined in our qudit paper
%or wiki page.
Motion_factor_state_raise_COM(:,:,1:2:qudit_level)=repmat(expm(-(-1^((qudit_level-1)/2)).*1i.*eta_COM.*(exp(1i.*w_COM.*t).*Motion_raise_COM+exp(-1i.*w_COM.*t).*Motion_lower_COM)),[1 1 numel(1:2:qudit_level)]);
Motion_factor_state_raise_COM(:,:,2:2:(qudit_level-1))=repmat(expm((-1^((qudit_level-1)/2)).*1i.*eta_COM.*(exp(1i.*w_COM.*t).*Motion_raise_COM+exp(-1i.*w_COM.*t).*Motion_lower_COM)),[1 1 numel(2:2:(qudit_level-1))]);
Motion_factor_state_lower_COM(:,:,1:2:qudit_level)=repmat(expm((-1^((qudit_level-1)/2)).*1i.*eta_COM.*(exp(1i.*w_COM.*t).*Motion_raise_COM+exp(-1i.*w_COM.*t).*Motion_lower_COM)),[1 1 numel(1:2:qudit_level)]);
Motion_factor_state_lower_COM(:,:,2:2:(qudit_level-1))=repmat(expm(-(-1^((qudit_level-1)/2)).*1i.*eta_COM.*(exp(1i.*w_COM.*t).*Motion_raise_COM+exp(-1i.*w_COM.*t).*Motion_lower_COM)),[1 1 numel(2:2:(qudit_level-1))]);
%Motion_factor_state_raise_Tilt(:,:,1:2:qudit_level)=repmat(expm(-(-1^((qudit_level-1)/2)).*1i.*eta_Tilt.*(exp(1i.*w_Tilt.*t).*Motion_raise_Tilt+exp(-1i.*w_Tilt.*t).*Motion_lower_Tilt)),[1 1 numel(1:2:qudit_level)]);
%Motion_factor_state_raise_Tilt(:,:,2:2:(qudit_level-1))=repmat(expm((-1^((qudit_level-1)/2)).*1i.*eta_Tilt.*(exp(1i.*w_Tilt.*t).*Motion_raise_Tilt+exp(-1i.*w_Tilt.*t).*Motion_lower_Tilt)),[1 1 numel(2:2:(qudit_level-1))]);
%Motion_factor_state_lower_Tilt(:,:,1:2:qudit_level)=repmat(expm((-1^((qudit_level-1)/2)).*1i.*eta_Tilt.*(exp(1i.*w_Tilt.*t).*Motion_raise_Tilt+exp(-1i.*w_Tilt.*t).*Motion_lower_Tilt)),[1 1 numel(1:2:qudit_level)]);
%Motion_factor_state_lower_Tilt(:,:,2:2:(qudit_level-1))=repmat(expm(-(-1^((qudit_level-1)/2)).*1i.*eta_Tilt.*(exp(1i.*w_Tilt.*t).*Motion_raise_Tilt+exp(-1i.*w_Tilt.*t).*Motion_lower_Tilt)),[1 1 numel(2:2:(qudit_level-1))]);
M_evol_array=nan(numel(psi),numel(psi),qudit_level-1);
for qq=1:qudit_level-1
    l=qq-ceil(qudit_level/2);
    M_evol_array(:,:,qq)=-1i.*Omega_bar.*cos(delta.*t).*...
        (kron(kron(((-1).^l).*1i.*State_raise(:,:,qq),State_I),kron(Motion_factor_state_raise_COM(:,:,qq),Motion_I_Tilt))...
        +kron(kron(State_I,((-1).^l).*1i.*State_raise(:,:,qq)),kron(Motion_factor_state_raise_COM(:,:,qq),Motion_I_Tilt))...
        +kron(kron(-((-1).^l).*1i.*State_lower(:,:,qq),State_I),kron(Motion_factor_state_lower_COM(:,:,qq),Motion_I_Tilt))...
        +kron(kron(State_I,-((-1).^l).*1i.*State_lower(:,:,qq)),kron(Motion_factor_state_lower_COM(:,:,qq),Motion_I_Tilt)));
end
M_evol=sum(M_evol_array,3);
dpsidt=M_evol*psi;