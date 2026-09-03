function subjectColorMap = getSubjectColorMap()
% getSubjectColorMap
%
% Returns a containers.Map that assigns a fixed, perceptually distinct
% RGB color to each subject ID.
%
% MATLAB: R2020b compatible
%
% Usage:
%   cmap = getSubjectColorMap();
%   rgb  = cmap('SP044');   % returns 1x3 RGB row vector

    % ---- Master subject list (fixed order = fixed colors) ----
    SUBJECT_LIST = { ...
        'SP035', ... %new
        'SP037', ... %new
        'SP043', ... %new
        'SP044', ...
        'SP046', ... %new
        'SP052', ... %new
        'SP053', ... %new
        'SP054', ...
        'SP058', ...
        'SP060', ...
        'SP061', ...
        'SP063', ...
        'SP065', ...
        'SP066', ...
        'SP067', ...
        'SP072', ...
        'SP075', ...
        'SP076'  ...
    };

    % ---- Base on tableau tab10 palette (RGB in 0..1) ----
    % Order: blue, (dark grey,) orange, green, red, purple,
    %        brown, pink, gray, olive, cyan
    tab11 = [ ...
        0.2000 0.2000 0.2000 %this is dark gray (SP035 is non canonical for now)
        0.2000 0.2000 0.2000 %this is dark gray (SP037 is non canonical for now)
        0.2000 0.2000 0.2000 %this is dark gray (SP043 is non canonical for now)
        0.1216 0.4667 0.7059
        0.2000 0.2000 0.2000 %this is dark gray (SP046 is non canonical for now)
        0.2000 0.2000 0.2000 %this is dark gray (SP052 is non canonical for now)
        0.2000 0.2000 0.2000 %this is dark gray (SP053 is non canonical for now)
        0.2000 0.2000 0.2000 %this is dark gray (SP054 is non canonical for now)
        1.0000 0.4980 0.0549
        0.1725 0.6275 0.1725
        0.8392 0.1529 0.1569
        0.5804 0.4039 0.7412
        0.2000 0.2000 0.2000 %this is dark gray (SP065 is non-learner)
        0.5490 0.3373 0.2941
        0.8902 0.4667 0.7608
        0.4980 0.4980 0.4980
        0.7373 0.7412 0.1333
        0.0902 0.7451 0.8118
    ];

    %add other options here if needed
    %...

    %make a choice here
    SUBJECT_COLORS = tab11;

    % ---- Build lookup map: subject ID -> RGB ----
    subjectColorMap = containers.Map( ...
        SUBJECT_LIST, ...
        num2cell(SUBJECT_COLORS, 2) ...
    );

end
