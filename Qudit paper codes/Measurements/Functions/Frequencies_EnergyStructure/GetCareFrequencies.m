function [Freqs] = GetCareFrequencies(G, LevelG, LevelP, CarrierFreq, GeomOrientation)
%Given a set of ground states occupied LevelG, and shelving states occupied
%LevelP, GetCareFrequencies calculates all of the frequencies that involve
%any of these
%Input G is the set of global variables, which contains atomic information
%G.I, G.J, G.Jp, G.EnergiesGround, G.EnergiesExcited, G.Fs, G.Fps, and
%G.MotionalFreq
%It outputs Freqs, which gives all transitions as rows of the form
%Col1: The actual frequency
%Col2:5: The F levels involved in this transition
%Col6: Clebschs Gordon Coefficient
%Col7: Order of the transition (0 is carrier)
%Col8: Whether or not it's an encoded transition, or a sideband of an
%encoded transition

%Retrieve all of the relevant frequencies given the orientation and energy
%structure
[Freqs, FreqsMot] = CalculateFreqs(G.I, G.J, G.Jp, CarrierFreq, ...
    G.EnergiesGround, G.EnergiesExcited, G.Fs, G.Fps, GeomOrientation, true, G.MotionalFreq);
%Only add MotionalFreqs to FreqsCare if it's nonempty and this error is
%turned on
if ~isempty(FreqsMot) && G.MotionalErrorOn
    Freqs = [Freqs; FreqsMot];
end

%Go through and get the encoded transitions indices and also get all of 
%the transitions involving encoded energy levels
EncodedList = [];
Care = [];
k = 1;
%Only go through and look at the levels
for i = 1:size(LevelG, 1)
    if GeomOrientation == "XZ" || GeomOrientation == "Average"
        %Ensure the list is always just an array, not a matrix
        EncodedList = [EncodedList reshape(find(all(Freqs(:, 2:5) == [LevelG(i, :) LevelP(i,:)], 2)).', 1, [])];
    elseif GeomOrientation == "Orthogonal"
        for j = 1:size(LevelP,1)
            if all(~isnan(LevelG(i, :))) && all(~isnan(LevelP(j, :))) && ...
                    abs(LevelG(i, 1) - LevelP(j, 1)) == 0 && abs(LevelG(i, 2) - LevelP(j, 2)) == 2
                EncodedList = [EncodedList reshape(find(all(Freqs(:, 2:5) == [LevelG(i, :) LevelP(j, :)], 2)).', 1, [])];
            end
        end
    end
    Care = [Care find(all(Freqs(:, 2:3) == LevelG(i, :), 2) | all(Freqs(:,4:5) == LevelP(i, :), 2)).'];
end
%Get rid of duplicates
Care = unique(Care);
[a, index] = intersect(Care, EncodedList);
Care(index) = [];
%Get frequencies that care about, adding column for whether encoded or not
FreqsCareEncoded = [Freqs(EncodedList, :) ones(length(EncodedList), 1)];
FreqsCareOther = [Freqs(Care, :) zeros(length(Care), 1)];
Freqs = [FreqsCareEncoded;FreqsCareOther];
end