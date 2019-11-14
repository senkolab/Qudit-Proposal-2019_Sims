g0_factor=0.5;
t_ends=thetas./(2*pi*Rabi_frequency);
t=0:deltat:sum(t_ends);
t_noise=0:0.5*deltat:sum(t_ends);
psi00=zeros(8,1);

uB=9.274009994*(10^(-24));
Planck_h=6.62607015*(10^(-34));

if qudit_level==5
    m_array=[-2 -1 0 1 2 -1 0 1];
    g_array=g0_factor.*[1 -1 1 -1 1 1 -1 1];
elseif qudit_level==3
    m_array=[-1 0 1 -2 -1 0 1 2];
    g_array=g0_factor.*[1 -1 1 1 -1 1 -1 1];
end

H_ideal=nan(numel(psi00),numel(psi00),1);
for pp=1:numel(thetas)
    [H_I,Off_resonant_freq]=Givens_Hamiltonian_TV(Rabi_frequency,spin_phases(pp),transitions(pp),g_array,m_array,off_resonant);
    if pp==1
        %H_I=repmat(H_I,[1 1 numel(t(t>=0 & t<=t_ends(pp)))]);
        H_I=repmat(H_I,[1 1 numel(t_noise(t_noise>=0 & t_noise<=t_ends(pp)))]);
    else
        %H_I=repmat(H_I,[1 1 numel(t(t>sum(t_ends(1:(pp-1))) & t<=sum(t_ends(1:(pp)))))]);
        H_I=repmat(H_I,[1 1 numel(t_noise(t_noise>sum(t_ends(1:(pp-1))) & t_noise<=sum(t_ends(1:(pp)))))]);
    end
    H_ideal=cat(3,H_ideal,H_I);
end
H_ideal(:,:,1)=[];

U_ideal=SingleQuditReconstruct(thetas,spin_phases,transitions,numel(psi00),0);

percent=0.1;
tic
error_array=nan(n_iter,1);
self_error_array=nan(n_iter,1);
for hh=1:n_iter
    psi00=zeros(8,1);
    %psi00(1:qudit_level)=[rand(qudit_level,1)+1i.*rand(qudit_level,1)];
    %psi00=psi00./sqrt(psi00'*psi00);
    psi00(1)=1;
    psi_ideal=U_ideal*psi00;
    psi=nan(numel(psi00),numel(t));
    psi(:,1)=psi00;
    %B_noise=Ornstein_Uhlenbeck(2000,2.7*(10^(-12)),t);
    B_noise=Ornstein_Uhlenbeck(2000,2.7*(10^(-12)),t_noise);
    B_noise=reshape(B_noise,[1 1 numel(B_noise)]);
    H_noise=2*pi.*repmat((10^(-6)).*(uB/Planck_h).*diag(g_array.*m_array),[1 1 numel(B_noise)])...
        .*repmat(B_noise,[numel(psi00) numel(psi00) 1]);
    H_net=H_ideal+H_noise;
for h=2:numel(t)
    %psi(:,h)=psi(:,h-1)-1i.*(H_net(:,:,h-1)*psi(:,h-1)).*(t(h)-t(h-1));
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