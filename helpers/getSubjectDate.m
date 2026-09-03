function [dt_bCW, dt_tCWFull] = getSubjectDate(subjectID)
% getSubjectDate
% Returns datetime of (1) first full task session (i.e. first biasedCW 
% session) and (2) first full contrast set session (i.e. corresponding to 
% training phase 5)

    [expMap1, expMap2] = getSubjectDateMap();

    if ~isKey(expMap1, subjectID)
        error('Unknown subject ID: %s', subjectID);
    end
    
    if ~isKey(expMap2, subjectID)
        error('Unknown subject ID: %s', subjectID);
    end

    dt_bCW = expMap1(subjectID);
    dt_tCWFull = expMap2(subjectID);

end