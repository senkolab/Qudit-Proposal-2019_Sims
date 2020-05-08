function Freq = CalculateFreqsHide(I, J, Energies, Fs, LevelG, LevelP)
%% Generate the F and mF levels from Fs, Fsp
F = [];
mF = [];
for i=1:length(Fs)
    F = [F ones(1, 2*Fs(i)+1)*Fs(i)];
    mF = [mF -Fs(i):1:Fs(i)];
end
%% Generate all of the frequency differences between all energy levels
FreqsLower = repmat(Energies, 1, length(Energies));
FreqsUpper = FreqsLower.';
FreqsTransitions = FreqsUpper - FreqsLower;
%Go through and pick out the relevant transitions from all of the
%transitions - these are those transitions that are allowed based on atomic
%rules and geometrical constraints from laser orientation and polarization
k = 0;
for i=1:length(mF)
    F0 = F(i);
    mF0 = mF(i);
    for j=1:length(mF)
        F1 = F(j);
        mF1 = mF(j);
        DeltaM = max(mF1, mF0) - min(mF1, mF0);
        if DeltaM > 1 || F1-F0 == 0
            continue
        end
        k = k + 1;
        RelevantFreqs(k,1:5) = [FreqsTransitions(i, j) F0 mF0 F1 mF1];
    end
end
RelevantFreqs = [abs(RelevantFreqs(:, 1)) RelevantFreqs(:, 2:5)];
%Order the matrix of relevant frequencies by absolute value of frequency
[temp, order] = sort(RelevantFreqs(:, 1));
Freqs = RelevantFreqs(order, :);
Care = find(all(Freqs(:, 2:3) == LevelG, 2) & all(Freqs(:,4:5) == LevelP, 2)).';
%Get rid of duplicates
Care = unique(Care);
Freq = Freqs(Care, :);
end