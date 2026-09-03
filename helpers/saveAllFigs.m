FolderName = 'C:\Users\Samuel\Documents\2PI\FENS2024\temp';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  FigNumber = get(FigHandle, 'Number');
  savefig(FigHandle, fullfile(FolderName, [num2str(FigNumber) '.fig']));
end