function [Probs, TransferTime] = Prob3level(Rabi, RabiMat, Linewidth, SweepMat, Sweep, FreqsCare3level, F, ClebschsCare3level, DecayTime)
    %ProbIdeal3Level is for tallying up the overall prob
    ProbIdeal3Level = ones(size(Rabi));
    %TotalTransferTime3 is for tallying up the overall transfer time
    TotalTransferTime3 = zeros(size(Rabi));
    %Go through each passage
    for i = 2:3
        %Get the prob of the this transfer
        prob3level = Prob5(Linewidth, RabiMat, SweepMat, FreqsCare3level, F, i, ClebschsCare3level);
        %Get the sweep rate with the best fidelity for each Rabi freq
        [ProbIdeal, index] = max(prob3level);
        Detunings = abs(FreqsCare3level(i) -FreqsCare3level);
        Detunings(i) = [];
        SmallestDetuning = min(Detunings);
        %Calculate the transfer time of this passage (optimal sweep rate; for each Rabi freq)
        TransferTime = 2*SmallestDetuning./Sweep(index);
        TransferTime = TransferTime.';
        %Add to transfer tally
        TotalTransferTime3 = TotalTransferTime3 + TransferTime;
        %Add in decay error
        ProbIdeal = ProbIdeal.*exp(-TransferTime/DecayTime);
        %Mult this transfer prob to the tally prob
        ProbIdeal3Level = ProbIdeal.*ProbIdeal3Level;
        %Add sweep data from this transfer to sweep data array
        %IdealSweeps3(i, :) = Sweep(index);
    end
    for i = 2:2
        %Get the prob of the this transfer
        prob3level = Prob5(Linewidth, RabiMat, SweepMat, FreqsCare3level, F, i, ClebschsCare3level);
        %Get the sweep rate with the best fidelity for each Rabi freq
        [ProbIdeal, index] = max(prob3level);
        Detunings = abs(FreqsCare3level(i) -FreqsCare3level);
        Detunings(i) = [];
        SmallestDetuning = min(Detunings);
        %Calculate the transfer time of this passage (optimal sweep rate; for each Rabi freq)
        TransferTime = 2*SmallestDetuning./Sweep(index);
        TransferTime = TransferTime.';
        %Add to transfer tally
        TransferTime = TotalTransferTime3 + TransferTime;
        %Add in decay error
        ProbIdeal = ProbIdeal.*exp(-TransferTime/DecayTime);
        %Mult this transfer prob to the tally prob
        Probs = ProbIdeal.*ProbIdeal3Level;
        %Add sweep data from this transfer to sweep data array
        %IdealSweeps3(i, :) = Sweep(index);
    end
end