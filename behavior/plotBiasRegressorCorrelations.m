function [Rcorr,regressors] = plotBiasRegressorCorrelations( ...
    trialsT,params,driftBias)
%PLOTBIASREGRESSORCORRELATIONS Compare full-fit bias regressors.
%
% INPUTS
%   trialsT   - IBL trials table (session or superSession)
%   params    - Struct with scalar fields:
%                 q            contrast exponent
%                 k            subjective block lag (trials)
%                 tauSide      stimulus-side history decay
%                 tauContrast  contrast-weighted history decay
%   driftBias - Optional Nx1 vector of full-fit drifting bias
%               (latent deviation from constant intercept).
%
% OUTPUTS
%   Rcorr      - Pearson correlation matrix
%   regressors - Nx3 or Nx4 matrix of raw bias regressors
%
% If driftBias is omitted, compares the three stimulus/block
% regressors only. Uses all trials with finite regressors.
%
% Lower triangle: scatter plots, linear fits, Pearson r
% Diagonal:       regressor distributions

if nargin<3
    driftBias = [];
end

n = height(trialsT);

%% Load stimulus and block variables

if ismember('contrastDiff',trialsT.Properties.VariableNames)
    X = double(trialsT.contrastDiff(:));
else
    cr = double(trialsT.contrastRight(:));
    cl = double(trialsT.contrastLeft(:));
    cr(isnan(cr)) = 0;
    cl(isnan(cl)) = 0;
    X = cr-cl;
end

pLeft = double(trialsT.probabilityLeft(:));

%% Identify session boundaries

if ismember('sessionID',trialsT.Properties.VariableNames)
    sessionCode = findgroups(trialsT.sessionID);
elseif ismember('sessionNr',trialsT.Properties.VariableNames)
    sessionCode = findgroups(trialsT.sessionNr);
else
    sessionCode = ones(n,1);
end

sessionCode = sessionCode(:);

starts = [1;find(diff(sessionCode)~=0)+1];
ends = [starts(2:end)-1;n];

%% 1. Delayed block regressor

B = zeros(n,1);
B(abs(pLeft-0.2)<1e-6) =  1;
B(abs(pLeft-0.8)<1e-6) = -1;

Bshift = zeros(n,1);

for s = 1:numel(starts)
    ids = (starts(s):ends(s))';

    idx = max(1,min(numel(ids), ...
        (1:numel(ids))'-params.k));

    Bshift(ids) = B(ids(idx));
end

%% 2-3. Exponentially weighted stimulus history

Hside = localHistory(X,params.q,params.tauSide, ...
    starts,ends,false);

Hcontrast = localHistory(X,params.q,params.tauContrast, ...
    starts,ends,true);

%% Assemble regressor matrix

regressors = [Bshift,Hside,Hcontrast];

names = {'Block','Stimulus side','Contrast weighted'};

colors = [
    0.00 0.45 0.70;   % Block
    0.00 0.62 0.45;   % Side
    0.85 0.37 0.01];  % Contrast

if ~isempty(driftBias)

    driftBias = double(driftBias(:));

    assert(numel(driftBias)==n, ...
        'driftBias must have one value per trial.');

    regressors = [regressors,driftBias];

    names{end+1} = 'Drifting bias';
    colors(end+1,:) = [0.58 0.40 0.74];
end

%% Select finite trials

keep = all(isfinite(regressors),2);

R = regressors(keep,:);

nReg = size(R,2);
nTrials = size(R,1);

assert(nTrials>=3,'Not enough valid trials.');

%% Correlations

Rcorr = corrcoef(R);

%% Subsample points for visualization only

maxScatter = 3000;

if nTrials>maxScatter
    scatterIdx = unique(round(linspace( ...
        1,nTrials,maxScatter)));
else
    scatterIdx = 1:nTrials;
end

%% Plot

figure('Color','w','Position',[100 100 900 850]);

tl = tiledlayout(nReg,nReg, ...
    'TileSpacing','compact','Padding','compact');

colPoints = [0.35 0.35 0.35];
colFit = [0.10 0.10 0.10];

for row = 1:nReg
    for col = 1:nReg

        ax = nexttile(tl,(row-1)*nReg+col);

        if row<col
            axis(ax,'off');
            continue
        end

        hold(ax,'on');

        x = R(:,col);
        y = R(:,row);

        if row==col

            % Diagonal: distribution
            histogram(ax,x,35, ...
                'FaceColor',colPoints, ...
                'EdgeColor','none', ...
                'FaceAlpha',0.65);

            %title(ax,names{row}, ...
            %    'Color',colors(row,:), ...
            %    'FontWeight','bold');

            if col==1
                ylabel(ax,'Count');
            end

        else

            % Scatter points (neutral grey)
            scatter(ax, ...
                x(scatterIdx),y(scatterIdx), ...
                8,colPoints, ...
                'filled', ...
                'MarkerFaceAlpha',0.12, ...
                'MarkerEdgeAlpha',0.12);

            % Linear fit over all trials
            if std(x)>0 && std(y)>0

                beta = polyfit(x,y,1);

                xFit = linspace(min(x),max(x),100);
                yFit = polyval(beta,xFit);

                plot(ax,xFit,yFit,'-', ...
                    'Color',colFit, ...
                    'LineWidth',1.5);

                r = Rcorr(row,col);

                text(ax,0.05,0.95, ...
                    sprintf('r = %.2f',r), ...
                    'Units','normalized', ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','top', ...
                    'FontSize',10, ...
                    'FontWeight','bold');
            end
        end

        %% Axis labels and formatting

        if row==nReg
            xlabel(ax,names{col}, ...
                'Color',colors(col,:), ...
                'FontWeight','bold');
        else
            ax.XTickLabel = [];
        end

        if row>col
            if col==1
                ylabel(ax,names{row}, ...
                    'Color',colors(row,:), ...
                    'FontWeight','bold');
            else
                ax.YTickLabel = [];
            end
        end

        box(ax,'off');
        ax.TickDir = 'out';
        ax.FontSize = 9;

    end
end

title(tl,sprintf( ...
    'Bias regressor correlations (full fit, N = %d trials)', ...
    nTrials));

end


%% Helper: normalized exponentially weighted stimulus history

function H = localHistory(X,q,tau,starts,ends,weighted)

n = numel(X);
H = zeros(n,1);

decay = exp(-1/tau);

for s = 1:numel(starts)

    numerator = 0;
    denominator = 0;

    for t = starts(s):ends(s)

        % History up to, but not including, current trial
        if denominator>0
            H(t) = numerator/denominator;
        end

        if ~isfinite(X(t)) || X(t)==0
            continue
        end

        if weighted
            w = abs(X(t))^q;
        else
            w = 1;
        end

        numerator = decay*numerator + w*sign(X(t));
        denominator = decay*denominator + w;
    end
end

end
