function [psychoM,RTs] = getPsycho_fitPseudosessions(trialsT,ks,nPseudo)
%GETPSYCHO_FITPSEUDOSESSIONS Psychometric fit with block null test.
%
% [psychoM,RTs] = getPsycho_fitPseudosessions(trialsT,ks,nPseudo)
%
% Fits the empirical PM model, then fits nPseudo randomized
% block sequences using identical fitting and lag selection.
%
% Adds to psychoM:
%   pseudoDevExplainedBlock : null distribution
%   pBlock                 : one-sided Monte Carlo p-value
%   significantBlock       : pBlock <= 0.05
%
% Defaults: ks=1:20, nPseudo=199

if nargin<2 || isempty(ks), ks = 1:20; end
if nargin<3 || isempty(nPseudo), nPseudo = 199; end

alpha = 0.05; %significance threshold

% Empirical model
[psychoM,RTs] = getPsycho_fit(trialsT,true,ks);

% Pseudo models
pseudoDE = nan(nPseudo,1);
pseudoBestLag = nan(nPseudo,1);
pseudoDB = nan(nPseudo,1);

for i = 1:nPseudo

    pseudoT = IBL_genPseudoBlocks(trialsT);

    pseudoM = getPsycho_fit(pseudoT,true,ks);

    pseudoDE(i) = pseudoM.devExplainedBlock;
    pseudoBestLag(i) = pseudoM.bestLag;
    pseudoDB(i) = pseudoM.deltabias;

end

% Require valid pseudo fits
if any(~isfinite(pseudoDE))
    error('Some pseudosession fits returned nonfinite deviance explained.');
end

% One-sided Monte Carlo significance test
pBlock = (1 + sum(pseudoDE >= psychoM.devExplainedBlock)) / (nPseudo+1);

% Store results
psychoM.pseudoDevExplainedBlock = pseudoDE;
psychoM.pseudoBestLag = pseudoBestLag;
psychoM.pseudoDeltabias = pseudoDB;

psychoM.pBlock = pBlock;
psychoM.significantBlock = pBlock <=alpha;

end