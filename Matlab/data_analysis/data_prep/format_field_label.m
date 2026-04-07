function label = format_field_label(name)
%FORMAT_FIELD_LABEL Convert field names (e.g., U3, S_11) to TeX-safe labels.
%   label = format_field_label('U3') -> 'U_{3}'
%   label = format_field_label('S_Mises') -> 'S_{\mathrm{Mises}}'

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

tok = regexp(name, '^([A-Za-z])_(\\d+)$', 'tokens', 'once');
if ~isempty(tok)
    label = sprintf('%s_{%s}', tok{1}, tok{2});
    return;
end

tok = regexp(name, '^([A-Za-z])_([A-Za-z]+)$', 'tokens', 'once');
if ~isempty(tok)
    if strcmpi(tok{2}, 'Mises')
        label = sprintf('%s_{\\mathrm{Mises}}', tok{1});
        return;
    end
    label = sprintf('%s_{%s}', tok{1}, tok{2});
    return;
end

label = strrep(name, '_', '\_');
label = strrep(label, '\', '\\');
end
