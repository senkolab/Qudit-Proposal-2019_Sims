function dpsidt=Qudit_Givens(t,psi,H,Off_resonant_freq,varargin)
%Input definition:
%t: Time.
%psi: Qudit state. 1D array.
%H: Rabi frequencies of Hamiltonian. 2D array.
%Off_resonant_freq: Time varying factors of Hamiltonian due to
%off-resonance.
%H_noise: Time-varying magnetic field noise component of the Hamiltonian.
%3D matrix.
%t_noise: Time array corresponding to H_noise. 1D array.

%Output definition:
%dpsidt: Rate of change of state. 1D array.
if numel(varargin)==2 %only used for time-varying Hamiltonian due to magnetic field noise
    H_noise=varargin{1};
    t_noise=varargin{2};
    H_noise_interp=interp3(1:1:size(H_noise,1),1:1:size(H_noise,2),t_noise,H_noise,1:1:size(H_noise,1),1:1:size(H_noise,2),t);
    Hnet=H+H_noise_interp;
    M_evol=-1i.*Hnet.*exp(1i.*Off_resonant_freq.*t);
    dpsidt=M_evol*psi;
elseif isempty(varargin)
    M_evol=-1i.*H.*exp(1i.*Off_resonant_freq.*t);
    dpsidt=M_evol*psi;
end