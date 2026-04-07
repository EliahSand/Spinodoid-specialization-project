function make_laminate(params)
%MAKE_LAMINATE Convenience entry point for laminate generation.
%   Uses the same parameter pattern as PSSCone: pass a struct with any
%   fields to override defaults in make_laminate_sheet.

if nargin < 1
    params = struct();
end

make_laminate_sheet(params);

end
