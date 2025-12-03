function label = format_field_label(name)
%FORMAT_FIELD_LABEL Convert field names (e.g., U3, S11) to TeX-safe labels.
%   label = format_field_label('U3') -> 'U_{3}'
%   label = format_field_label('SMises') -> 'S_{\mathrm{Mises}}'

arguments
    name (1, :) char
end

if isempty(name)
    label = '';
    return;
end

prefix = name(1);
suffix = name(2:end);
if ~isempty(suffix) && all(isstrprop(suffix, 'digit'))
    if isstrprop(prefix, 'alpha')
        label = sprintf('%s_{%s}', prefix, suffix);
        return;
    end
end

label = strrep(name, '_', '\_');
label = strrep(label, '\', '\\');
end
