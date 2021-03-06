function Care = GetCareTransitions(Level, FreqsAbs, GeomOrientation)
%This function finds which transitions we care about (meaning they involve
%encoded levels) given a particular encoding. 
%The output is in the form of a vector specifying which indices of the input
%frequency matrix contains the transitions we care about
%The inputs are: 
%Level: How many levels are being encoded?
%FreqsAbs: a list of all transitions in the form of 1:6:Freq, Fs, Clebschs.
%GeomOrientation: what Geometric orientation we're using (XZ, Orthogonal,
%Average)
%Hide: allows us to hide an encoded level elsewhere, because
%of motional sidebands or other reasons. This is a vector where the first
%entry tells which encoded level to hide, and the next tell 
G = getGlobals();
%Encodings for 3- 5- and 7-level qudits most times
Encoded3n = [...
    1 -1 1 -1;...
    1 1 1 1;...
    2 0 2 0];
Encoded5n = [...
    Encoded3n;...
    2 -2 2 -2;...
    2 2 2 2];
Encoded7n = [...
    Encoded3n;...
    1 0 1 0;...
    2 1 2 1;...
    2 -1 2 -1;...
    2 2 2 2];

EncodedBase3 = [...
    1 -1;...
    1 1;...
    2 0];
EncodedBase5 = [...
    EncodedBase3;...
    1 -1;...
    1 1;...
    2 0];
EncodedBase7 = [...
    EncodedBase3;...
    1 0;...
    2 1;...
    2 -1;...
    2 2];
if GeomOrientation == "XZ" || GeomOrientation == "Average"
    EncodedBaseE3 = EncodedBase3;
    EncodedBaseE5 = EncodedBase5;
    EncodedBaseE7 = EncodedBaseE7;
elseif GeomOrientation == "Orthogonal"
    EncodedBaseE3 = [...
        1 1;...
        1 -1;...
        2 2];
    EncodedBaseE5 = [...
        EncodedBaseE3;...
        2 0;...
        2 0];
    EncodedBaseE7 = [...
        EncodedBaseE3;...
        1 0;...
        2 -1;...
        2 1;...
        2 0];
end

if GeomOrientation == "XZ" || GeomOrientation == "Average"
    Encoded3 = Encoded3n;
    Encoded5 = Encoded5n;
    Encoded7 = Encoded7n;
elseif GeomOrientation == "Orthogonal"
    Encoded3 = [...
        1 -1 1 1;...
        1 1 1 -1;...
        2 0 2 2];
    Encoded5 = [...
        Encoded3n;...
        2 -2 2 0;...
        2 2 2 0];
    Encoded7 = [...
        Encoded3n;...
        1 0 1 0;...
        2 1 2 -1;...
        2 -1 2 1;...
        2 2 2 0];
end

%Encodings for 3- 5- and 7-level qudits 
if Level == 3
    EncodedMatrix = Encoded3;
elseif Level == 5
    EncodedMatrix = Encoded5;
elseif Level == 7
    EncodedMatrix = Encoded7;
end
%Go through and get the encoded transitions indices and also get all of 
%the transitions involving encoded energy levels
Encoded = [];
Care = [];
k = 1;
for i = 1:size(EncodedMatrix, 1)
    Encoded = [Encoded find(all(FreqsAbs(:, 2:5) == EncodedMatrix(i, :), 2) & FreqsAbs(:, 7) == 0)];
    Care = [Care find(all(FreqsAbs(:, 2:3) == EncodedMatrix(i, 1:2), 2) | all(FreqsAbs(:,4:5) == EncodedMatrix(i, 3:4), 2)).'];
    k = k + 1;
end
%Get rid of duplicates
Care = unique(Care);
[a, index] = intersect(Care, Encoded);
Care(index) = [];
%Add the encoded levels at the beginning
Care = [Encoded Care];
end