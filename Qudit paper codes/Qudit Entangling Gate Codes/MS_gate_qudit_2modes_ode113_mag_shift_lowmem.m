%Omega_bar as defined in this code is twice the value of Omega_bar as defined in the paper
%due to a different way of generalizing the MS gate when this code was initially written.

tic %start timer for code

%Initialize qudit levels and MS gate parameters
%qudit_level=5;
s=(qudit_level-1)/2;
theta=pi/4; %Geometric phase shift for MS gate
w_COM=2.*pi.*fq_COM; %Motional frequency
w_Tilt=2.*pi.*fq_Tilt;
%delta=2.*pi.*5.01; %Laser detuning
delta=w_COM+n_cycle*2*pi/gate_time;
%eta_COM=0.1; %Lamb-Dicke parameter
eta_tilt=eta_COM*sqrt(w_COM/w_Tilt);

%Generalized S_x operator for qudit
S_operator_array=zeros(qudit_level,qudit_level,qudit_level);
for h=1:(qudit_level-1)
    l=h-ceil(qudit_level/2);
    S_operator_array(h+1,h,h)=sqrt(s.*(s+1)-l.*(l+1))./2;
    S_operator_array(h,h+1,h)=sqrt(s.*(s+1)-l.*(l+1))./2;
end
S_operator=sum(S_operator_array,3);
[eigVec,eigVal]=eig(S_operator); %Get eigenvalues and eigenvectors of generalized S_x operator

%Initialize input qudit states
Input_qudit1=zeros(qudit_level,1);
Input_qudit2=zeros(qudit_level,1);
Input_qudit1(end)=1; %Input qudit 1
Input_qudit2(end)=1; %Input qudit 2
Input_state=kron(Input_qudit1,Input_qudit2);

Input_qudits_x=kron(eigVec',eigVec')*Input_state; %Input qudits in x basis

%Compute phase shifts for each eigenvectors
phase_array=zeros(numel(Input_qudits_x),1);
for h=1:qudit_level
    phase_array((qudit_level*(h-1)+1):(qudit_level*h))=sign(w_COM-delta).*theta.*(eigVal(h,h)+diag(eigVal)).^2;
end
%Compute output state of MS gate after gaining geometric phases
Output_qudits_x=exp(1i.*phase_array).*Input_qudits_x;
Output_qudits=kron(eigVec,eigVec)*Output_qudits_x;

%Initialize S+ and S- operators for qudits
State_raise=zeros(qudit_level,qudit_level,qudit_level);
State_lower=zeros(qudit_level,qudit_level,qudit_level);
for h=1:(qudit_level-1)
    l=h-ceil(qudit_level/2);
    State_raise(h+1,h,h)=sqrt(s.*(s+1)-l.*(l+1))./2;
    State_lower(h,h+1,h)=sqrt(s.*(s+1)-l.*(l+1))./2;
end

%Initialize motional operators
N_number_COM=20;
Motion_I_COM=diag(ones(N_number_COM,1));
Motion_raise_COM=diag(sqrt(1:1:N_number_COM-1),-1);
Motion_lower_COM=diag(sqrt(1:1:N_number_COM-1),1);

N_number_Tilt=2;
Motion_I_Tilt=diag(ones(N_number_Tilt,1));
Motion_raise_Tilt=diag(sqrt(1:1:N_number_Tilt-1),-1);
Motion_lower_Tilt=diag(sqrt(1:1:N_number_Tilt-1),1);

%Initialize motional state
n_COM=0:1:N_number_COM-1;
n_COM=n_COM';
%n_average_COM=0.1;
%Pn_COM=exp(-n_average_COM).*(n_average_COM.^n_COM)./factorial(n_COM);
Pn_COM=(n_average_COM.^n_COM)./((n_average_COM+1).^(n_COM+1));
%Motion_state_COM=sqrt(Pn_COM);

n_Tilt=0:1:N_number_Tilt-1;
n_Tilt=n_Tilt';
n_average_Tilt=0;
%Pn_Tilt=exp(-n_average_Tilt).*(n_average_Tilt.^n_Tilt)./factorial(n_Tilt);
Pn_Tilt=(n_average_Tilt.^n_Tilt)./((n_average_Tilt+1).^(n_Tilt+1));
%Motion_state_Tilt=sqrt(Pn_Tilt);

Ln=nan(size(n_COM));
Ln_1=nan(size(n_COM));
for h=1:numel(Ln)
    m=0:1:n_COM(h);
    Ln(h)=sum(((-1).^m).*(factorial(n_COM(h)+1)./(factorial(n_COM(h)-m).*factorial(1+m))).*(eta_COM.^(2.*m))./factorial(m));
end
Ln_1(2:end)=Ln(1:end-1);
Ln_1_sqn=[0; (Ln_1(2:end).^2)./n_COM(2:end)];
Gn=1-exp(-eta_COM.^2).*((Ln.^2)./(n_COM+1)-Ln_1_sqn);
%Gn=0;
Omega_bar=(1/sqrt(n_cycle))*sqrt(1+sum(Pn_COM.*Gn.*(1-Gn))./sum(Pn_COM.*(1-Gn).^2)).*sqrt((w_COM+delta).*theta.*((w_COM-delta).^2)./(pi.*w_COM.*(eta_COM.^2))); %Rabi frequency

%Input state with both qudit states and motional state
%psi=kron(Input_state,kron(Motion_state_COM,Motion_state_Tilt));

%Initialize time array and time steps for numerical simulation
%{
deltat=0.001;
t_space=0.01;
t_space=deltat*round(t_space/deltat);
t_max=gate_time;
t_max=t_space*round(t_max/t_space);
t_array=linspace(0,t_max,round(t_max/t_space)+1);
t=0;
psi_t=nan(numel(psi),numel(t_array));
psi_t(:,1)=psi;
%}
State_I=diag(ones(qudit_level,1));

%Initialize arrays for numerical simulation
%{
Motion_factor_state_raise_COM=nan(N_number_COM,N_number_COM,qudit_level);
Motion_factor_state_lower_COM=nan(N_number_COM,N_number_COM,qudit_level);
Motion_factor_state_raise_Tilt=nan(N_number_Tilt,N_number_Tilt,qudit_level);
Motion_factor_state_lower_Tilt=nan(N_number_Tilt,N_number_Tilt,qudit_level);
M_evol_array=nan((qudit_level^2)*N_number_COM*N_number_Tilt,(qudit_level^2)*N_number_COM*N_number_Tilt,qudit_level);
%}

if qudit_level==3
    State_mag_shift = 2.*pi.*diag([37.74039 -0.124146 -37.74039]./2);
    State_mag_shift = State_mag_shift./(10^9);
elseif qudit_level==5
    State_mag_shift = 2.*pi.*diag([-75.66723 37.74039 -0.124146 -37.74039 75.66723]./2);
    State_mag_shift = State_mag_shift./(10^9);
else
    sprintf('Undefined qudit levels');
end

if mag_shift_on==0
    State_mag_shift=0.*State_mag_shift;
end

disp('Initialization complete')
fidel=0;

for c=1:N_number_COM
    for cc=1:N_number_Tilt
        disp(['n_COM = ' num2str(c-1) ', n_Tilt = ' num2str(cc-1)])
        if Pn_COM(c)*Pn_Tilt(cc)<=0.00001
            continue
        end
        Motion_state_COM=(n_COM==c-1);
        Motion_state_Tilt=(n_Tilt==cc-1);
        psi=kron(Input_state,kron(Motion_state_COM,Motion_state_Tilt));
        
        if heattime~=0
            [tm,psi_t]=ode113(@(t,psi) MS_gate_fun_mag_shift(t,psi,qudit_level,gate_time,n_cycle,eta_COM,w_COM,w_Tilt,Omega_bar,...
                State_raise,State_lower,State_I,Motion_raise_COM,Motion_lower_COM,...
                Motion_raise_Tilt,Motion_lower_Tilt,State_mag_shift),[0 heattime*gate_time],psi,options);
            tm1=tm;
            psi_t1=psi_t;
            psi=psi_t(end,:);
            psi=psi.';
            psi=kron(kron(State_I,State_I),kron(diag(ones(N_number_COM-1,1),-1),Motion_I_Tilt))*psi;
        end
        
        [tm,psi_t]=ode113(@(t,psi) MS_gate_fun_mag_shift(t,psi,qudit_level,gate_time,n_cycle,eta_COM,w_COM,w_Tilt,Omega_bar,...
                State_raise,State_lower,State_I,Motion_raise_COM,Motion_lower_COM,...
                Motion_raise_Tilt,Motion_lower_Tilt,State_mag_shift),[heattime*gate_time gate_time],psi,options);
        
        tm2=tm;
        psi_t2=psi_t;
        
        %Compute fidelity
        %fidel=fidel=sum(abs(kron(Output_qudits',kron(Motion_I_COM,Motion_I_Tilt))*psi_t(end,:).').^2);;
        %fidel=sum(abs(kron(Output_qudits',kron(Motion_I_COM,Motion_I_Tilt))*U_tm(:,:,end)*psi).^2);
        psi_end=psi_t(end,:).';
        %psi_end=U_tm(:,:,end)*psi;
        rho=psi_end*psi_end';
        %rho=U_tm(:,:,end)*rho_in*U_tm(:,:,end)';
        rho_state_COM=nan((qudit_level^2)*N_number_COM,(qudit_level^2)*N_number_COM);
        for h=1:size(rho_state_COM,1)
            for hh=1:size(rho_state_COM,1)
                rho_state_COM(h,hh)=sum(diag(rho((N_number_Tilt*(h-1)+1):h*N_number_Tilt,(N_number_Tilt*(hh-1)+1):hh*N_number_Tilt)));
            end
        end
        rho_state=nan(qudit_level^2,qudit_level^2);
        for h=1:size(rho_state,1)
            for hh=1:size(rho_state,1)
                rho_state(h,hh)=sum(diag(rho_state_COM((N_number_COM*(h-1)+1):h*N_number_COM,(N_number_COM*(hh-1)+1):hh*N_number_COM)));
            end
        end
        fidel2=Output_qudits'*rho_state*Output_qudits;
        fidel=fidel+Pn_COM(c).*Pn_Tilt(cc).*fidel2;
    end
end
toc