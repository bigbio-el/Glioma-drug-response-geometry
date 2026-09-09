function inspect_glioma_loom
[fileName,filePath]=uigetfile({'*.loom','Loom files (*.loom)'},'Select glioma loom file');
if isequal(fileName,0), return; end
f=fullfile(filePath,fileName);
fprintf('\nFILE: %s\n\n',f);
mi=h5info(f,'/matrix'); disp(mi.Dataspace.Size);
for groupName = ["/row_attrs","/col_attrs"]
    fprintf('\n%s\n',groupName);
    info=h5info(f,char(groupName));
    for i=1:numel(info.Datasets)
        p=char(groupName + "/" + string(info.Datasets(i).Name));
        fprintf('  %s',p);
        try
            x=h5read(f,p);
            sz=size(x); fprintf('   size='); fprintf('%dx',sz); fprintf('\b ');
            if iscell(x)
                q=string(x(:));
                fprintf(' examples: %s',strjoin(q(1:min(4,end)), ' | '));
            elseif isstring(x)
                q=x(:); fprintf(' examples: %s',strjoin(q(1:min(4,end)), ' | '));
            elseif ischar(x)
                fprintf(' char');
            else
                fprintf(' %s',class(x));
            end
        catch ME
            fprintf(' [read error: %s]',ME.message);
        end
        fprintf('\n');
    end
end
end
