function T = compactROItable(AllROIs, opts)
%COMPACTALLROIS Compact a large 1xN struct array of ROI records into a lean table.
%
%   T = compactROItable(AllROIs)
%   T = compactROItable(AllROIs, opts)
%
% Key behavior:
%   - Drops tstat_pseudo
%   - Converts to table
%   - Tightens numeric types
%   - Converts pos -> Nx3 single array
%   - Keeps subject/date/session as separate fields

    arguments
        AllROIs (1,:) struct
        opts.dropTstatPseudo (1,1) logical = true
        opts.categoricalSSD (1,1) logical = true
        opts.dateToDatetime (1,1) logical = false
        opts.clearInput (1,1) logical = false
    end

    % ---- 1) Drop heavy field
    if opts.dropTstatPseudo && isfield(AllROIs, 'tstat_pseudo')
        AllROIs = rmfield(AllROIs, 'tstat_pseudo');
    end

    % ---- 2) Convert to table
    T = struct2table(AllROIs);

    % ---- 3) subject / date / session handling
    if ismember('subject', T.Properties.VariableNames) && opts.categoricalSSD
        T.subject = categorical(T.subject);
    end
    if ismember('session', T.Properties.VariableNames) && opts.categoricalSSD
        T.session = categorical(T.session);
    end
    if ismember('date', T.Properties.VariableNames)
        if opts.dateToDatetime
            if ~isdatetime(T.date)
                T.date = datetime(T.date,'InputFormat','yyyy-MM-dd');
            end
        elseif opts.categoricalSSD
            T.date = categorical(T.date);
        end
    end

    % ---- 4) Tighten numeric storage
    castIfExists('iCell',  'uint32');
    castIfExists('iROI',   'uint32');
    castIfExists('iFOV',   'uint16');
    castIfExists('annot',  'uint32');
    castIfExists('ccf_id', 'uint32');

    castIfExists('tstat_empirical','single');
    castIfExists('p',               'single');

    if ismember('h', T.Properties.VariableNames)
        T.h = logical(T.h);
    end

    % ---- 5) pos -> Nx3 single array
    if ismember('pos', T.Properties.VariableNames)
        try
            P = vertcat(T.pos);   % N x 3
            if size(P,2) ~= 3
                error('Expected pos to be Nx3 after vertcat.');
            end
            T.pos = single(P);
        catch ME
            warning('compactAllROIs:posConvert', ...
                'Could not convert pos to Nx3 single array: %s', ME.message);
        end
    end

    % ---- 6) Optional: clear input in caller workspace
    if opts.clearInput
        evalin('caller','clear AllROIs');
    end

    % Helper
    function castIfExists(varName, typeName)
        if ismember(varName, T.Properties.VariableNames) ...
                && ~isa(T.(varName), typeName)
            T.(varName) = feval(typeName, T.(varName));
        end
    end
end
