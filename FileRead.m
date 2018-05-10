function [ InpTextInrow ] = FileRead( InpName,UpLimRowNum )

  %InpRead: To read all context in the file by row, and a cell will be returned.
  % -----------------------------------------------------
  % Input:
  % InpName: A String,the file name to be read.
  % UpLimRowNum: A integer, equal to or more than the largest row number.
  % Output:
  % InpTextInrow: A cell, containing all the text in the file.
  % -----------------------------------------------------

  Fpn = fopen (InpName, 'rt');   % Open the files.
  ii = 0; % Initialize a index for the row number.
  InpTextInrow = cell(UpLimRowNum,1);  % Initialize a cell to store the text read from the file.
  while feof(Fpn) ~= 1  % When the pointer comes to the end of the file, Fpn = 1; otherwise, Fpn = 0.
    ii = ii + 1;          % All rows in the file will be read.
    InpTextInrow{ii,1} = fgetl(Fpn);  % Read the next row.
  end
  if ii < UpLimRowNum
    InpTextInrow = InpTextInrow(1:ii,1);  %Remove empty row in InptextInrow
  end
  fclose(Fpn);

end