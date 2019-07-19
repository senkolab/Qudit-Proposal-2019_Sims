function U=SingleQuditReconstruct(thetas,spin_phases,transitions,totallevels,theta_G)
%Input definition:
%thetas: Gate phases of the pulses in radian. 1D array.
%spin_phases: Spin phases of the pulses in radian. 1D array.
%transitions: Transitions of the pulses. 1D array of integers.
%totallevels: Total states in a qudit state, including unencoded state.
%Integer
%theta_G: Global phase in radian.

%Output definition:
%U: Time evolution matrix reconstructed from a sequence of pulses as
%specified in input. 2D matrix.
U=diag(ones(totallevels,1));
for h=1:numel(thetas)
    U_dummy=diag(ones(size(U,1),1));
    U_dummy(transitions(h),transitions(h))=cos(thetas(h)/2);
    U_dummy(transitions(h)+1,transitions(h)+1)=cos(thetas(h)/2);
    U_dummy(transitions(h)+1,transitions(h))=-1i*sin(thetas(h)/2)*exp(1i*spin_phases(h));
    U_dummy(transitions(h),transitions(h)+1)=-1i*sin(thetas(h)/2)*exp(-1i*spin_phases(h));
    U=U_dummy*U;
end
U=exp(1i*theta_G).*U;