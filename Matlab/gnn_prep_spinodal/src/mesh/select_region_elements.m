function region = select_region_elements(inpData, varargin)
%SELECT_REGION_ELEMENTS Select elements/nodes for target region from an Abaqus mesh.
%
% region = SELECT_REGION_ELEMENTS(inpData, ...)
%
% Name-Value options:
%   'ElsetName'         : explicit elset name to use (default '')
%   'AutoDetectElset'   : auto-pick spinodal/top elset if explicit is empty (default true)
%   'PatternPriority'   : ordered tokens used for auto-detect
%                         default {'spinodal','top'}
%
% Output region fields:
%   source               - 'elset_explicit' | 'elset_auto' | 'all_elements'
%   elset_name           - selected elset name or ''
%   element_labels       - selected element labels
%   element_connectivity - selected element node labels (cell)
%   node_labels          - unique node labels in selected elements

    p = inputParser;
    p.addParameter('ElsetName', '', @(x) ischar(x) || isstring(x));
    p.addParameter('AutoDetectElset', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PatternPriority', {'spinodal', 'top'}, @(x) iscell(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;

    requestedName = char(string(opts.ElsetName));
    autoDetect = logical(opts.AutoDetectElset);
    patterns = cellstr(string(opts.PatternPriority));

    selectedElementLabels = inpData.element_labels;
    source = 'all_elements';
    chosenElset = '';

    if ~isempty(inpData.elset_names)
        idx = [];

        if ~isempty(strtrim(requestedName))
            idx = find(strcmpi(inpData.elset_names, requestedName), 1);
            if isempty(idx)
                error('select_region_elements:ElsetNotFound', ...
                    'Requested elset "%s" not found in inp.', requestedName);
            end
            source = 'elset_explicit';
            chosenElset = inpData.elset_names{idx};
        elseif autoDetect
            idx = autodetect_elset(inpData.elset_names, patterns);
            if ~isempty(idx)
                source = 'elset_auto';
                chosenElset = inpData.elset_names{idx};
            end
        end

        if ~isempty(idx)
            selectedElementLabels = unique(inpData.elset_values{idx}(:), 'stable');
        end
    end

    [isMember, loc] = ismember(selectedElementLabels, inpData.element_labels);
    loc = loc(isMember & loc > 0);
    loc = unique(loc, 'stable');

    selectedElementLabels = inpData.element_labels(loc);
    selectedConnectivity = inpData.element_connectivity(loc);

    selectedNodeLabels = zeros(0, 1);
    if ~isempty(selectedConnectivity)
        selectedNodeLabels = unique(vertcat(selectedConnectivity{:}), 'stable');
    end

    region = struct();
    region.source = source;
    region.elset_name = chosenElset;
    region.element_labels = selectedElementLabels;
    region.element_connectivity = selectedConnectivity;
    region.node_labels = selectedNodeLabels;
end

function idx = autodetect_elset(elsetNames, patterns)
    lowerNames = lower(string(elsetNames));
    idx = [];

    for p = 1:numel(patterns)
        pat = lower(string(patterns{p}));
        matches = find(contains(lowerNames, pat));
        if ~isempty(matches)
            idx = matches(1);
            return;
        end
    end
end

