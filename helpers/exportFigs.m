function exportFigs(outDir, varargin)
% exportFigs
%
% Usage:
%   export_all_figures_svg(outDir)
%   export_all_figures_svg(outDir, 'overwrite', true)
%   export_all_figures_svg(outDir, 'overwrite', true, 'mysuffix') 
%
% Rules:
%   - If last argument is char/string → treated as suffix
%   - Otherwise no suffix
%
% Options:
%   'overwrite'  (false)
%   'timestamp'  (false)

% ---- extract suffix (must be last) ----
suffix = '';

if ~isempty(varargin)
    lastArg = varargin{end};

    if ischar(lastArg) || isstring(lastArg)
        suffix = char(string(lastArg));
        varargin = varargin(1:end-1); % remove suffix
    end
end

% ---- parse name-value pairs ----
p = inputParser;
p.addParameter('overwrite', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('timestamp', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

opt.overwrite = logical(opt.overwrite);
opt.timestamp = logical(opt.timestamp);

% ---- sanitize suffix ----
if ~isempty(suffix)
    suffix = regexprep(suffix, '[^\w\d-_ ]', '');
    suffix = strrep(suffix, ' ', '_');
end

% ---- ensure output dir ----
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figs = findall(0, 'Type', 'figure');

if isempty(figs)
    warning('No open figures to export.');
    return
end

for i = 1:numel(figs)
    fig = figs(i);

    % base name
    nm = get(fig, 'Name');
    if isempty(nm)
        nm = sprintf('Figure_%d', fig.Number);
    end

    nm = regexprep(nm, '[^\w\d-_ ]', '');
    nm = strrep(nm, ' ', '_');

    % append suffix
    if ~isempty(suffix)
        nm = [nm '_' suffix];
    end

    fname = fullfile(outDir, [nm '.svg']);

    % overwrite handling
    if isfile(fname)
        if opt.timestamp
            ts = datestr(now, 'yyyymmdd_HHMMSS');
            fname = fullfile(outDir, [nm '_' ts '.svg']);
        elseif ~opt.overwrite
            warning('File exists, skipping: %s', fname);
            continue
        end
    end

    % export
    try
        exportgraphics(fig, fname, 'ContentType', 'vector');
        %exportgraphics(fig, fname, 'ContentType', 'image', 'Resolution', 300);
        fprintf('Saved: %s\n', fname);
    catch ME
        warning('Failed to save figure "%s": %s', nm, ME.message);
    end
end
end