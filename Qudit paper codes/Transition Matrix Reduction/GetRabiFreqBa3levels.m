P12_F=[1 2]; %F states of P_1/2 level
P32_F=[0 1 2 3]; %F states of P_3/2 level

State0=[2 -1]; %[F m_F] levels of state 0
State1=[1 0]; %[F m_F] levels of state 1
State2=[2 1]; %[F m_F] levels of state 2

T1=[0 0 0 0]; %Initialize array which contains transition matrix element 
%reduction coefficients for the transition 0 to 1. First and second
%elements correspond to coupling to P_1/2 level with different polazation
%combinations. Third and fourth elements correspond to coupling to P_3/2
%level with different polarization combinations.
T2=[0 0 0 0]; %Initialize array which contains transition matrix element 
%reduction coefficients for the transition 1 to 2. First and second
%elements correspond to coupling to P_1/2 level with different polazation
%combinations. Third and fourth elements correspond to coupling to P_3/2
%level with different polarization combinations.

J=[1/2 1/2 3/2 3/2]; %J levels of P states
q1=[0 1 0 1]; %Polarization of first beam
q2=[-1 0 -1 0]; %Polarization of second beam

%Compute the reduction coefficients for 0 to 1 transition
for h=1:numel(T1)
    if h==1 || h==2
        P_F=P12_F;
    else
        P_F=P32_F;
    end
    for hh=1:numel(P_F)
        for hhh=1:(2*P_F(hh)+1)
            mfp=-P_F(hh)-1+hhh;
            T1(h)=T1(h)+FmtoLCoefficient(P_F(hh),mfp,2,-1,1,q1(h),J(h),1/2,3/2,1,0,1/2)*FmtoLCoefficient(P_F(hh),mfp,1,0,1,q2(h),J(h),1/2,3/2,1,0,1/2);
        end
    end
end

%Compute the reduction coefficients for 1 to 2 transition
for h=1:numel(T2)
    if h==1 || h==2
        P_F=P12_F;
    else
        P_F=P32_F;
    end
    for hh=1:numel(P_F)
        for hhh=1:(2*P_F(hh)+1)
            mfp=-P_F(hh)-1+hhh;
            T2(h)=T2(h)+FmtoLCoefficient(P_F(hh),mfp,1,0,1,q1(h),J(h),1/2,3/2,1,0,1/2)*FmtoLCoefficient(P_F(hh),mfp,2,1,1,q2(h),J(h),1/2,3/2,1,0,1/2);
        end
    end
end