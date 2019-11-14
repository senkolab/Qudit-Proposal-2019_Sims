function U=SingleQuditReconstruct(thetas,spin_phases,transitions,totallevels,theta_G)
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