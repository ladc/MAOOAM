function namelist=readin_nml(varargin)
% This function reads any !standard! fortran namelist file. the output is a
% structure with fieldnames corresponding to the namelist variables.
for i=1:length(varargin)
    fid = fopen(varargin{i});
    if( fid == -1 )
        error(['Unable to open ' fname]);
    end
    while( 1 )
        tline = fgetl(fid);
        if( isequal(tline,-1) )
            break;
        end
        try
            k=strfind(tline, '=');
            e=strfind(tline, '!');
            if isempty(e); e=length(tline)+1;end;
            var_name=strtrim(tline(1:k-1));
            if any(strfind(var_name,':'));
                eval(['namelist.',var_name,'=[',strtrim(tline(k+1:e-1)),'];']);
            else
                eval(['namelist.',var_name,'=',strtrim(tline(k+1:e-1)),';']);
            end
        catch
        end
    end
    fclose(fid);
end