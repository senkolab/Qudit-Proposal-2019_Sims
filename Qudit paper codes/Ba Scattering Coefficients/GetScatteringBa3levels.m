P12_F=[1 2]; %F levels of the P_1/2 level
P32_F=[0 1 2 3]; % F levels of the P_3/2 level

State0=[2 -1]; %[F m_F] levels of state 0
State1=[1 0]; %[F m_F] levels of state 1
State2=[2 1]; %[F m_F] levels of state 2
States=[State0; State1; State2]; %Array of the states

%Initialize transition matrix element reduction coefficient arrays for all
%the polarizations
rplus=zeros(size(States,1),2);
rminus=zeros(size(States,1),2);
b0=zeros(size(States,1),2);

J=[1/2 3/2]; %J levels of P_1/2 and P_3/2 states

%Compute the transition matrix element reduction coefficients for
%sigma-plus polarization, first column of rplus corresponds to transition
%to P_1/2 state while second column corresponds to P_3/2. Each row
%corresponds to State 0 to State 2 in increasing order. Each element in the
%array is summed over the the transition to all F and m_F levels in a P
%state.
for h1=1:size(States,1); 
    for h2=1:numel(J)
        if h2==1
            P_F=P12_F;
        else
            P_F=P32_F;
        end
        for h4=1:numel(P_F)
            for h5=1:(2*P_F(h4)+1)
                mfp=-P_F(h4)-1+h5;
                rplus(h1,h2)=rplus(h1,h2)+FmtoLCoefficient(P_F(h4),mfp,States(h1,1),States(h1,2),1,1,J(h2),1/2,3/2,1,0,1/2)^2;
            end
        end
    end
end

%Compute the transition matrix element reduction coefficients for
%sigma-minus polarization, first column of rplus corresponds to transition
%to P_1/2 state while second column corresponds to P_3/2. Each row
%corresponds to State 0 to State 2 in increasing order. Each element in the
%array is summed over the the transition to all F and m_F levels in a P
%state.
for h1=1:size(States,1);
    for h2=1:numel(J)
        if h2==1
            P_F=P12_F;
        else
            P_F=P32_F;
        end
        for h4=1:numel(P_F)
            for h5=1:(2*P_F(h4)+1)
                mfp=-P_F(h4)-1+h5;
                rminus(h1,h2)=rminus(h1,h2)+FmtoLCoefficient(P_F(h4),mfp,States(h1,1),States(h1,2),1,-1,J(h2),1/2,3/2,1,0,1/2)^2;
            end
        end
    end
end

%Compute the transition matrix element reduction coefficients for
%pi polarization, first column of rplus corresponds to transition
%to P_1/2 state while second column corresponds to P_3/2. Each row
%corresponds to State 0 to State 2 in increasing order. Each element in the
%array is summed over the the transition to all F and m_F levels in a P
%state.
for h1=1:size(States,1);
    for h2=1:numel(J)
        if h2==1
            P_F=P12_F;
        else
            P_F=P32_F;
        end
        for h4=1:numel(P_F)
            for h5=1:(2*P_F(h4)+1)
                mfp=-P_F(h4)-1+h5;
                b0(h1,h2)=b0(h1,h2)+FmtoLCoefficient(P_F(h4),mfp,States(h1,1),States(h1,2),1,0,J(h2),1/2,3/2,1,0,1/2)^2;
            end
        end
    end
end