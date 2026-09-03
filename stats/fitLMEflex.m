function lme = fitLMEflex(tbl, dv, ivs, groupVar, varargin)
% dv:        char/string, dependent variable name in tbl
% ivs:       cellstr/string array of predictor names
% groupVar:  char/string, random-effects grouping variable in tbl
% opts:      struct with fields (all optional):
%   .RandomSlopeVars  = cellstr of IVs to include as random slopes (default: {})
%   .IncludeIntercept = true/false (default: true)
%   .Categoricals     = cellstr of variable names to cast to categorical (default: {})
%   .FitMethod        = 'REML' (default) or 'ML'
%
% Example:
%   lme = fitLMEflex(tbl,'AUC',{'DeltaBias','Age'},'Subject', ...
%                    struct('RandomSlopeVars',{{'DeltaBias'}}, 'Categoricals',{{'Subject'}}));

    % ---- Defaults
    opts = struct( ...
        'RandomSlopeVars', {{}}, ...
        'IncludeIntercept', true, ...
        'Categoricals',    {{}}, ...
        'FitMethod',       'REML');

    % ---- Parse varargin (struct OR name-value)
    if ~isempty(varargin)
        if isstruct(varargin{1})
            s = varargin{1};
            f = fieldnames(s);
            for k = 1:numel(f)
                opts.(f{k}) = s.(f{k});
            end
        else
            % name–value pairs
            for k = 1:2:numel(varargin)
                name = varargin{k};
                val  = varargin{k+1};
                opts.(name) = val;
            end
        end
    end

    % ensure categoricals if requested
    for c = opts.Categoricals
        vn = string(c);
        if ~iscategorical(tbl.(vn))
            tbl.(vn) = categorical(tbl.(vn));
        end
    end

    % fixed part 
    if isempty(ivs)
        if opts.IncludeIntercept
            fixed = "1";
        else
            fixed = "-1";
        end
    else
        fixedIVs = strjoin(ivs, ' + ');
        if opts.IncludeIntercept
            fixed = "1 + " + fixedIVs;
        else
            fixed = "-1 + " + fixedIVs;
        end
    end

    % random part
    if isempty(opts.RandomSlopeVars)
        randPart = sprintf('(1|%s)', groupVar);
    else
        randIVs = strjoin(["1", string(opts.RandomSlopeVars)], ' + ');
        randPart = sprintf('(%s|%s)', randIVs, groupVar);
    end

    % full formula
    formula = sprintf('%s ~ %s + %s', dv, fixed, randPart);

    % fit
    lme = fitlme(tbl, formula, 'FitMethod', opts.FitMethod);
end
