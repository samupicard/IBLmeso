function centile = comp_percentile(distribution, value)
%COMP_PERCENTILE Conservative rank-based percentile.
%   percentile = comp_percentile(distribution, value)
%
%   Computes the percentile of "value" relative to a finite sample
%   "distribution", including the value itself in the sample.
%
%   Uses midrank tie handling and the conservative formula:
%       p = (r_mid - 0.5) / n
%
%   Output is in percent (0–100), but never exactly 0 or 100.

    % Ensure column vector
    distribution = distribution(:);

    % Include the value in the sample
    S = [distribution; value];
    n = numel(S);

    % Count strictly less than
    L = sum(S < value);

    % Count equal (includes the appended value itself)
    E = sum(S == value);

    % Midrank
    r_mid = L + (E + 1) / 2;

    % Conservative percentile (strictly between 0 and 1)
    p = (r_mid - 0.5) / n;

    % Convert to percent
    centile = 100 * p;

end