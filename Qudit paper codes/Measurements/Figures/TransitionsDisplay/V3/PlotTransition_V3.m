function [Plot, Line, Text, Area] = PlotTransition_V3(Trans, Color, MaxClebschs, Detuning, GeomOrientation)
%This function plots a single transition in a line with a marker
%The outputs are each of the graphical elements created by this plot,
%including the marker (in plot), the line, the text denoting the transition
%F numbers, and the area denoting the laser sweep
G = getGlobals_V3();
Freq = Trans(1);
Height = 1;
if G.SameHeight
    Height = Height*G.SameHeightHeight;
else
    %Scale the height of this line by its Clebschs Gordan Coefficient
    Height = Height*abs(Trans(6))/MaxClebschs;
end
%Carrier transition?
if Trans(7) == 0
    MarkerColor = Color;
    ShowText = true;
    %Encoded transition?
    if Trans(8) == 1
        ShowArea = true;
        LineWidth = G.EncodedLineWidth;
        LineStyle = G.EncodedLineStyle;
        Alpha = G.EncodedAlpha;
%         if GeomOrientation == "XZ" || GeomOrientation == "Average"
%             %TextString = ['F = ' num2str(Trans(2)) ', m_F = ' num2str(Trans(3))];
%             TextString = sprintf("(%i, %i)\leftrightarrow (%i, %i)'", num2str(Trans(2)), num2str(Trans(3)), num2str(Trans(4)), num2str(Trans(5)));
%         elseif GeomOrientation == "Orthogonal"
%             TextString = ['F = ' num2str(Trans(2)) ', m_F = ' num2str(Trans(3)) ', F'' = ' num2str(Trans(4)) ', m_{F''} = ' num2str(Trans(5))];
%         end
    %Other transition?
    else
        ShowArea = false;
        %TextString = ['F = ' num2str(Trans(2)) ', m_F = ' num2str(Trans(3)) ', F'' = ' num2str(Trans(4)) ', m_{F''} = ' num2str(Trans(5))];
        LineWidth = G.OtherLineWidth;
        LineStyle = G.OtherLineStyle;
        Alpha = G.OtherAlpha;
    end
    TextString = sprintf('(%+i, %+i)\\leftrightarrow (%+i, %+i)''', Trans(2), Trans(3), Trans(4), Trans(5));
    %Motional sideband transition?
else
    %We don't need to show the area or text for motional sidebands
    ShowArea = false;
    ShowText = false;
    %1st order motional sideband?
    if abs(Trans(7)) == 1
        LineWidth = G.MotionalOneLineWidth;
        LineStyle = G.MotionalOneLineStyle;
        Alpha = G.MotionalOneAlpha;
        MarkerColor = G.MotionalOneMarkerColor;
        Height = Height*G.MotionalHeightOne;
    %2nd order motional sideband?
    elseif abs(Trans(7)) == 2
        LineWidth = G.MotionalTwoLineWidth;
        LineStyle = G.MotionalTwoLineStyle;
        Alpha = G.MotionalTwoAlpha;
        MarkerColor = G.MotionalOneMarkerColor;
        Height = Height*G.MotionalHeightTwo;
    end
end
%Plot the line at the desired frequency
Line = line([Freq Freq], [0 Height], ...
    'Color', Color, ...
    'LineWidth', LineWidth, ...
    'LineStyle', LineStyle);
%Plot the marker at the end of this line
Plot = plot([Freq Freq], [Height Height], '^', ...
    'MarkerEdge', Color, ...
    'MarkerFace', MarkerColor);
%Show the transition F numbers
if G.ShowRelevantText && ShowText
    Text = text(Freq+G.TextOffSet, G.TextHeight+Height, ...
        TextString, ...
        'Color', Color, ...
        'FontSize', G.TextFontSize);
    set(Text,'Rotation',90);
else
    Text = false;
end
%Show the area the laser sweeps out
if ShowArea
    Area = area([Freq-Detuning Freq+Detuning], ...
        [Height*G.SweepHeight Height*G.SweepHeight], ...
        'FaceColor', Color, ...
        'LineStyle', 'none', ...
        'FaceAlpha', G.SweepAlpha);
else
    Area = false;
end
end