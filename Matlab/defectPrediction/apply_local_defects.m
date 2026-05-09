function [spinLayer, applied] = apply_local_defects(spinLayer, Svox, L, defectSpec)
%APPLY_LOCAL_DEFECTS Stamp cracks and/or holes into the top spin layer.
%
%   spinLayer : logical [N x N x tsV], [row=y, col=x, z] convention (post-permute).
%   Svox      : voxel size in metres (x and y; square grid assumed).
%   L         : physical side length of the periodic tile (metres).
%   defectSpec: scalar struct or struct array, each entry with fields:
%                 type       — 'crack' or 'hole'
%                 position   — [x_m, y_m] physical centre
%                 (crack)  length_m, theta_deg
%                 (hole)   radius_m, count

N    = size(spinLayer, 1);
tsV  = size(spinLayer, 3);
nD   = numel(defectSpec);

applied = repmat(struct('type','','ix0',0,'iy0',0,'extra',[]), nD, 1);

% Square-spiral offsets (unit vectors) for multi-hole placement.
spiralUnits = [1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; 1 -1; -1 -1];

for di = 1:nD
    d    = defectSpec(di);
    x_m  = d.position(1);
    y_m  = d.position(2);
    ix0  = mod(round(x_m / Svox), N);   % 0-based for modular arithmetic
    iy0  = mod(round(y_m / Svox), N);

    applied(di).type = d.type;
    applied(di).ix0  = ix0 + 1;
    applied(di).iy0  = iy0 + 1;

    if strcmp(d.type, 'crack')
        theta  = deg2rad(d.theta_deg);
        nSteps = ceil(d.length_m / Svox);
        half   = (nSteps - 1) / 2;

        if isfield(d, 'width_m') && d.width_m > Svox
            half_w = floor(round(d.width_m / Svox) / 2);
        else
            half_w = 0;
        end

        for k = 0:(nSteps - 1)
            t  = k - half;
            dx = round(t * cos(theta));
            dy = round(t * sin(theta));
            for pw = -half_w:half_w
                px = round(pw * (-sin(theta)));
                py = round(pw *   cos(theta));
                ix = mod(ix0 + dx + px, N) + 1;
                iy = mod(iy0 + dy + py, N) + 1;
                spinLayer(iy, ix, :) = false;
            end
        end

    elseif strcmp(d.type, 'hole')
        r_vox  = d.radius_m / Svox;
        count  = d.count;
        centres = zeros(count, 2);   % 0-based
        centres(1,:) = [ix0, iy0];
        for ci = 2:count
            off = spiralUnits(ci - 1, :);
            centres(ci,:) = [ix0 + round(off(1) * 2.5 * r_vox), ...
                             iy0 + round(off(2) * 2.5 * r_vox)];
        end

        bR = ceil(r_vox);
        for ci = 1:count
            cx = centres(ci, 1);
            cy = centres(ci, 2);
            for dy = -bR:bR
                for dx = -bR:bR
                    if dx^2 + dy^2 <= r_vox^2
                        ix = mod(cx + dx, N) + 1;
                        iy = mod(cy + dy, N) + 1;
                        spinLayer(iy, ix, :) = false;
                    end
                end
            end
        end
    end
end

% Enforce periodic-boundary tileability after all stamps.
spinLayer(end,:,:) = spinLayer(1,:,:);
spinLayer(:,end,:) = spinLayer(:,1,:);
end
