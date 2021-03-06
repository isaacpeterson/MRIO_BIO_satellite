function rv=binread(fn, precision)
%
% Written by Keiichiro Kanemoto and Dan Moran
% Last revised: April 3 2013.
%

% read count data
[path, name, ext]=fileparts(fn);

% if exist([addsep(path) 'count_' name '.txt'], 'file'),
% 	sz=dlmread([addsep(path) 'count_' name '.txt'],'\t');
% else
%  % DM changed this for octave compatibility on 2012 5 30
%    if exist('OCTAVE_VERSION'),
% 		if ~strfind(fn,'sparse') & ~strfind(fn,'full'), warning(['No count_ file found for ' name]);
% 		else  sz=dlmread([addsep(path) 'count_' name '.txt'],'\t'); end;
% 		
%    else
% 		if strcmpi('.bin', fn(end-3:end)),  error(['No count_ file found for ' name]); end; 
%     end
% %    if ~strfind(fn,'sparse') & ~strfind(fn,'full'), warning(['No count_ file found for ' name]); end;
% 	sz = Inf;
% end;
% 
% set precision

if nargin<2, precision = 'double'; end;
% 
% read binary data
if ~exist(fn,'file'), error(['binread() says: No such file ' fn]);  end;
    
fid = fopen(fn, 'r');
% 
if strcmpi(ext,'.bin'),
    if length(sz) <= 2, %vector
        rv = fread(fid,sz,precision,0,'ieee-le');
    else %matrix
        rv = fread(fid,Inf,precision,0,'ieee-le');
        rv = reshape(rv, sz);
    end
elseif strcmpi(ext,'.spbin')
    ndim = fread(fid,1,precision,0,'ieee-le');
    issp = fread(fid,1,precision,0,'ieee-le');
    precis = fread(fid,1,precision,0,'ieee-le'); % precision. 8, 4, or 1 bytes. Currently always 8 (double)
    if ndim == 1, %vector
        if issp,
            szdim1 = fread(fid,1,precision,0,'ieee-le'); % get the sz of dimension 1, so we can osrt out whether to retun a row or column vector
            len = fread(fid,1,precision,0,'ieee-le');
            lennz = fread(fid,1,precision,0,'ieee-le');
            nzinds = fread(fid,lennz,'uint32',0,'ieee-le');
            nzvals = fread(fid,lennz,'double',0,'ieee-le');
            rv = zeros(len,1);
            rv(nzinds) = nzvals;
            if szdim1==1, rv = reshape(rv, 1,len); end;
        else
            szdim1 = fread(fid,1,precision,0,'ieee-le'); % get the sz of dimension 1, so we can osrt out whether to retun a row or column vector
            len = fread(fid,1,precision,0,'ieee-le');
            rv = fread(fid,len,precision,0,'ieee-le');
            if szdim1==1, rv = reshape(rv, 1,len); end;
        end;
        
    elseif ndim == 2, %matrix
        sz = fread(fid,2,'double',0,'ieee-le');
        if issp,
            spvec = fread(fid,Inf,precision,0,'ieee-le');
            spvec = reshape(spvec, length(spvec)/3,3);
            rv = sptensor(double(spvec(:, 1:2)), double(spvec(:, 3)), double(max(spvec(:, 1:2))));
            
        else
            
            rv = fread(fid,Inf,precision,0,'ieee-le');
            rv = reshape(rv, sz');
        end;
            
    else % tensor
        error([fn ' is corrupt. First byte should be 1 or 2, depending on whether data is vector or matrix']);
    end
   
else
    error(['binread can only read .bin/.spbin files. Not: ' name ext]);
end
fclose(fid);
end
