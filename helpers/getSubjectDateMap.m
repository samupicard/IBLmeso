function [FullTaskDateMap, FullContrastSetDateMap] = getSubjectDateMap()
% getSubjectDateMap
%
% Returns a set of containers.Map mapping subject ID to (1) first full task
% session (i.e. first biasedCW session) and (2) first full contrast set
% session (i.e. corresponding to training phase 5)
%
% Usage:
%   expMap               = getSubjectDateMap();
%   dt_bCW               = expMap('SP044');
%
%   or
%   [expMap1, expMap2]   = getSubjectDateMap();
%   dt_bCW               = expMap1('SP044');  
%   dt_tCWFull           = expMap2('SP044');  



    % ---- Master subject list (must match your global subject list) ----
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

    % ---- First biasedChoiceWorld session dates ----
    FullTask_dates_str = { ...
        '2023-02-15', ... % SP035 ... %new
        '2023-02-14', ... % SP037 ... %new
        '2023-06-16', ... % SP043 ... %new
        '2023-08-02', ... % SP044
        '2023-11-24', ... % SP046 ... %new
        '2024-05-14', ... % SP052 ... %new
        '2024-02-15', ... % SP053 ... %new
        '2024-02-16', ... % SP054
        '2024-07-12', ... % SP058
        '2024-07-29', ... % SP060
        '2024-12-10', ... % SP061
        '2025-03-17', ... % SP063
        '2024-12-01', ... % SP065 (non-learner, this is dummy date)
        '2025-04-08', ... % SP066
        '2025-06-03', ... % SP067
        '2025-08-21', ... % SP072
        '2025-12-03', ... % SP075
        '2025-11-10'  ... % SP076
    };

    % ---- First full contrast set session dates ----
    FullContrastSet_dates_str = { ...
        '2023-01-01', ... % SP035 %TODO
        '2023-01-01', ... % SP037 %TODO
        '2023-01-01', ... % SP043 %TODO
        '2023-01-01', ... % SP044 %TODO
        '2023-01-01', ... % SP046 %TODO
        '2024-01-01', ... % SP052 %TODO
        '2024-01-01', ... % SP053 %TODO
        '2024-01-01', ... % SP054 %TODO
        '2024-01-01', ... % SP058 %TODO
        '2024-07-15', ... % SP060
        '2024-12-10', ... % SP061 %TODO
        '2025-01-01', ... % SP063 %TODO
        '2024-12-01', ... % SP065 (non-learner, this is dummy date)
        '2025-01-01', ... % SP066 %TODO
        '2025-01-01', ... % SP067 %TODO
        '2025-01-01', ... % SP072 %TODO
        '2025-01-01', ... % SP075 %TODO
        '2025-01-01'  ... % SP076 %TODO
    };

    % ---- Parse to datetime ----
    dateFmt = 'yyyy-MM-dd';
    FullTask_dates = datetime(FullTask_dates_str, 'InputFormat', dateFmt);
    FullContrastSet_dates = datetime(FullContrastSet_dates_str, 'InputFormat', dateFmt);

    % ---- Build lookup maps ----
    FullTaskDateMap = containers.Map(SUBJECT_LIST, num2cell(FullTask_dates));
    FullContrastSetDateMap = containers.Map(SUBJECT_LIST, num2cell(FullContrastSet_dates));


end
