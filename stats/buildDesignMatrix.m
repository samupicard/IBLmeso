function X = buildDesignMatrix(x1, x2, x3)
    % Inputs are column vectors
    x1x3 = x1 .* x3;
    x2x3 = x2 .* x3;

    % Add a column of ones for the intercept
    X = [ones(length(x1), 1), x1, x2, x3, x1x3, x2x3];
end
