function inp = read_abaqus_inp(inpPath)
%READ_ABAQUS_INP Parse Abaqus .inp mesh file (nodes, elements, elsets).
%
% inp = READ_ABAQUS_INP(inpPath)
%
% Returns a struct with fields:
%   path                 - source .inp path
%   node_labels          - Nx1 Abaqus node labels
%   node_coords          - Nx3 node coordinates
%   element_labels       - Mx1 element labels
%   element_connectivity - Mx1 cell, each cell is element node labels
%   element_types        - Mx1 cell with element type (e.g., 'S4R')
%   elset_names          - Kx1 cell of elset names
%   elset_values         - Kx1 cell of element labels in each elset
%
% Notes:
% - Supports *Node, *Element, *Elset blocks.
% - Supports *Elset, generate syntax.
% - Ignores comment lines starting with '**'.

    arguments
        inpPath (1, :) char
    end

    if ~isfile(inpPath)
        error('read_abaqus_inp:FileNotFound', 'File not found: %s', inpPath);
    end

    fid = fopen(inpPath, 'r');
    if fid < 0
        error('read_abaqus_inp:OpenFailed', 'Could not open file: %s', inpPath);
    end
    c = onCleanup(@() fclose(fid));

    nodeLabels = zeros(0, 1);
    nodeCoords = zeros(0, 3);

    elementLabels = zeros(0, 1);
    elementConnectivity = cell(0, 1);
    elementTypes = cell(0, 1);

    elsetNameToIndex = containers.Map('KeyType', 'char', 'ValueType', 'double');
    elsetNames = cell(0, 1);
    elsetValues = cell(0, 1);

    mode = '';
    currentElementType = '';
    currentElsetName = '';
    currentElsetGenerate = false;
    lineNumber = 0;

    while true
        rawLine = fgetl(fid);
        if ~ischar(rawLine)
            break;
        end
        lineNumber = lineNumber + 1;

        line = strtrim(rawLine);
        if isempty(line) || startsWith(line, '**')
            continue;
        end

        if startsWith(line, '*')
            [keyword, attrs] = parse_keyword(line);
            switch keyword
                case 'NODE'
                    mode = 'NODE';

                case 'ELEMENT'
                    mode = 'ELEMENT';
                    if isfield(attrs, 'TYPE')
                        currentElementType = attrs.TYPE;
                    else
                        currentElementType = '';
                    end

                case 'ELSET'
                    mode = 'ELSET';
                    if ~isfield(attrs, 'ELSET')
                        error('read_abaqus_inp:BadElset', ...
                            'Missing elset name at line %d in %s', lineNumber, inpPath);
                    end
                    currentElsetName = attrs.ELSET;
                    currentElsetGenerate = isfield(attrs, 'GENERATE');

                    if ~isKey(elsetNameToIndex, currentElsetName)
                        elsetNames{end + 1, 1} = currentElsetName; %#ok<AGROW>
                        elsetValues{end + 1, 1} = zeros(0, 1); %#ok<AGROW>
                        elsetNameToIndex(currentElsetName) = numel(elsetNames);
                    end

                otherwise
                    mode = '';
            end
            continue;
        end

        switch mode
            case 'NODE'
                vals = parse_numeric_csv_line(line, false);
                if numel(vals) < 2
                    continue;
                end
                nodeLabel = round(vals(1));
                xyz = zeros(1, 3);
                nCoord = min(3, max(0, numel(vals) - 1));
                if nCoord > 0
                    xyz(1:nCoord) = vals(2:1 + nCoord);
                end

                nodeLabels(end + 1, 1) = nodeLabel; %#ok<AGROW>
                nodeCoords(end + 1, :) = xyz; %#ok<AGROW>

            case 'ELEMENT'
                vals = parse_numeric_csv_line(line, false);
                if isempty(vals)
                    continue;
                end

                % Abaqus may continue long element connectivity on lines that start with comma.
                isContinuation = startsWith(strtrim(rawLine), ',');
                if isContinuation
                    if isempty(elementConnectivity)
                        error('read_abaqus_inp:BadElementContinuation', ...
                            'Element continuation before first element at line %d in %s', ...
                            lineNumber, inpPath);
                    end
                    elementConnectivity{end} = [elementConnectivity{end}; round(vals(:))];
                else
                    if numel(vals) < 2
                        continue;
                    end
                    elementLabels(end + 1, 1) = round(vals(1)); %#ok<AGROW>
                    elementConnectivity{end + 1, 1} = round(vals(2:end).'); %#ok<AGROW>
                    elementTypes{end + 1, 1} = currentElementType; %#ok<AGROW>
                end

            case 'ELSET'
                vals = round(parse_numeric_csv_line(line, true));
                if isempty(vals)
                    continue;
                end
                idx = elsetNameToIndex(currentElsetName);
                if currentElsetGenerate
                    expanded = expand_generate(vals, lineNumber, inpPath, currentElsetName);
                    elsetValues{idx} = [elsetValues{idx}; expanded(:)];
                else
                    elsetValues{idx} = [elsetValues{idx}; vals(:)];
                end
        end
    end

    for i = 1:numel(elsetValues)
        if ~isempty(elsetValues{i})
            elsetValues{i} = unique(elsetValues{i}, 'stable');
        end
    end

    inp = struct();
    inp.path = inpPath;
    inp.node_labels = nodeLabels;
    inp.node_coords = nodeCoords;
    inp.element_labels = elementLabels;
    inp.element_connectivity = elementConnectivity;
    inp.element_types = elementTypes;
    inp.elset_names = elsetNames;
    inp.elset_values = elsetValues;
end

function [keyword, attrs] = parse_keyword(line)
    body = strtrim(line(2:end));
    parts = strsplit(body, ',');
    keyword = upper(strtrim(parts{1}));
    attrs = struct();

    if numel(parts) <= 1
        return;
    end

    for i = 2:numel(parts)
        token = strtrim(parts{i});
        if isempty(token)
            continue;
        end

        eqPos = strfind(token, '=');
        if isempty(eqPos)
            key = matlab.lang.makeValidName(upper(token));
            attrs.(key) = true;
            continue;
        end

        splitPos = eqPos(1);
        key = matlab.lang.makeValidName(upper(strtrim(token(1:splitPos - 1))));
        value = strtrim(token(splitPos + 1:end));
        attrs.(key) = value;
    end
end

function vals = parse_numeric_csv_line(line, integersOnly)
    if nargin < 2
        integersOnly = false;
    end

    parts = strsplit(line, ',');
    parts = parts(~cellfun('isempty', strtrim(parts)));
    if isempty(parts)
        vals = zeros(0, 1);
        return;
    end

    vals = str2double(parts(:));
    vals = vals(~isnan(vals));
    if integersOnly
        vals = round(vals);
    end
end

function expanded = expand_generate(vals, lineNumber, inpPath, elsetName)
    if mod(numel(vals), 3) ~= 0
        error('read_abaqus_inp:BadGenerate', ...
            ['Invalid *Elset, generate line at line %d in %s for elset %s. ', ...
             'Expected multiples of 3 values (start, stop, step).'], ...
            lineNumber, inpPath, elsetName);
    end

    expanded = zeros(0, 1);
    for i = 1:3:numel(vals)
        a = vals(i);
        b = vals(i + 1);
        s = vals(i + 2);

        if s == 0
            error('read_abaqus_inp:BadGenerateStep', ...
                'Zero step in *Elset, generate at line %d in %s', lineNumber, inpPath);
        end

        if (b >= a && s < 0) || (b <= a && s > 0)
            s = -s;
        end

        expanded = [expanded; (a:s:b).']; %#ok<AGROW>
    end
end

