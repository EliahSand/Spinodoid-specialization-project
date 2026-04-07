function label = format_ratio_label(tag)
%FORMAT_RATIO_LABEL Convert trXX labels to r = value strings.
%   format_ratio_label('tr200') -> 'r = 2'

arguments
    tag (1, :) char
end

digits = regexp(tag, '\d+', 'match', 'once');
if isempty(digits)
    label = tag;
    return;
end

ratio = str2double(digits) / 100;
if abs(round(ratio) - ratio) < 1e-9
    label = sprintf('r = %d', round(ratio));
else
    label = sprintf('r = %.2g', ratio);
end
end
