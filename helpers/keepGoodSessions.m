function Good = keepGoodSessions(All, goodsessions)

    fn = fieldnames(All);
    nSess = numel(All.sessID);   % or size of some session-dependent field

    for i = 1:numel(fn)
        name = fn{i};
        val  = All.(name);

        % skip non-array stuff (like strings, etc.)
        if ~isnumeric(val) && ~islogical(val) && ~iscell(val)
            Good.(name) = val;
            continue;
        end

        % if last dimension is "session", index it
        if size(val, ndims(val)) == nSess
            idx = repmat({':'}, 1, ndims(val));
            idx{end} = goodsessions;      % logical indexing on last dim
            Good.(name) = val(idx{:});
        else
            % keep as-is (e.g. PM_fitXs, global params, etc.)
            Good.(name) = val;
        end
    end
end
