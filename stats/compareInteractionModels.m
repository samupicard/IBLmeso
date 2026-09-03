function [F_vals, p_vals, R2_improvement] = compareInteractionModels(x1, x2, x3, Y)
% x1, x2, x3: column vectors [N x 1]
% Y: matrix of multiple responses [N x M]
% Returns:
%   F_vals: [1 x M] F-statistics
%   p_vals: [1 x M] p-values for each comparison
%   R2_improvement: [1 x M] improvement in R² from adding interactions

    [n, m] = size(Y);

    % Create design matrices
    X_base = [ones(n,1), x1, x2, x3];                 % 4 predictors
    X_full = [X_base, x1.*x3, x2.*x3];                % 6 predictors

    % Compute coefficients
    B_base = X_base \ Y;
    B_full = X_full \ Y;

    % Predictions
    Y_hat_base = X_base * B_base;
    Y_hat_full = X_full * B_full;

    % Residuals and RSS
    RSS_base = sum((Y - Y_hat_base).^2, 1);  % [1 x M]
    RSS_full = sum((Y - Y_hat_full).^2, 1);  % [1 x M]

    % Degrees of freedom
    q = size(X_full, 2) - size(X_base, 2);  % q = 2
    df_full = n - size(X_full, 2);          % df = n - 6

    % Compute F-statistic
    numerator = (RSS_base - RSS_full) / q;
    denominator = RSS_full / df_full;
    F_vals = numerator ./ denominator;

    % Compute p-values using F-distribution
    p_vals = 1 - fcdf(F_vals, q, df_full);

    % Optional: R² improvement
    TSS = sum((Y - mean(Y)).^2, 1);
    R2_base = 1 - RSS_base ./ TSS;
    R2_full = 1 - RSS_full ./ TSS;
    R2_improvement = R2_full - R2_base;
end