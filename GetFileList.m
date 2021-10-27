function [files,nf]=GetFileList(dirn,stype)
% [files,nf]=GetFileList(dirn,stype);
% INPUT (strings in single-quotes)
%   dirn is a directory name, such as 'C:\Active\Heinsohn\'  % Must end '\'
%    to read in the current directory, use '' (empty)
%   stype is the kind of file, such as '.ttt' or '.Master.transmission'
% OUTPUT  
%   files is a list of the files (without the directory in front)
%   nf    is the number of files (can be 0)
% Note: will recognize .Master.tra with stype='.tra'
drp=[dirn '*' stype];
drl=dir(drp); [nf,f]=size(drl);
files=[];
for f=1:nf
  temp=drl(f).name;
  files=[files; cellstr(temp)];
end;
end %get file list