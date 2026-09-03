function myCustomPlot(varargin)
%MYCUSTOMPLOT Plot data with optional parent axes and title, using positional args.
%
%   myCustomPlot(data)
%   myCustomPlot(parent, data)
%   myCustomPlot(parent, data, titleStr)

    p = inputParser;
    p.FunctionName = 'myCustomPlot';

    % Optional parent argument (can be omitted)
    addOptional(p, 'parent', [], @(x) isempty(x) || isgraphics(x, {'axes', 'uipanel', 'tiledlayout'}));

    % Required data argument (expected after parent)
    addOptional(p, 'data', [], @(x) isstruct(x) && all(isfield(x, {'x', 'y1', 'y2'})));

    % Optional title
    addOptional(p, 'titleStr', '', @(x) ischar(x) || isstring(x));

    % Parse inputs
    parse(p, varargin{:});
    r = p.Results;

    % Create axes if parent is not provided
    if isempty(r.parent)
        figure;
        ax1 = axes('Position', [0.1 0.6 0.8 0.3]);
        ax2 = axes('Position', [0.1 0.1 0.8 0.3]);
    else
        ax1 = axes('Parent', r.parent, 'Position', [0.1 0.6 0.8 0.3]);
        ax2 = axes('Parent', r.parent, 'Position', [0.1 0.1 0.8 0.3]);
    end

    % Plot
    plot(ax1, r.data.x, r.data.y1, 'b');
    plot(ax2, r.data.x, r.data.y2, 'r');
    linkaxes([ax1, ax2], 'x');

    % Add title if given
    if ~isempty(r.titleStr)
        sgtitle(r.titleStr);
    end
end
