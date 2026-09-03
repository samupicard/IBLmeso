function r = columnwisePearson(A, B)
A = A - mean(A, 1, 'omitnan');
B = B - mean(B, 1, 'omitnan');

denom = sqrt(sum(A.^2, 1, 'omitnan') .* sum(B.^2, 1, 'omitnan'));
r = sum(A .* B, 1, 'omitnan') ./ denom;

r(denom == 0) = NaN;
end