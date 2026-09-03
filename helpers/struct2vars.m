function struct2vars(s)
    f = fieldnames(s);
    for k = 1:numel(f)
        assignin('caller', f{k}, s.(f{k}));
    end
end