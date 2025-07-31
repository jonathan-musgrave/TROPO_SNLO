
function write_allStructs_to_ascii(allStructs, filename)
    % Writes all numeric fields from structs in allStructs to an ASCII text file

    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    % Loop over each struct in allStructs
    structNames = fieldnames(allStructs);
    for i = 1:numel(structNames)
        sName = structNames{i};
        s = allStructs.(sName);
        fieldNames = fieldnames(s);

        for j = 1:numel(fieldNames)
            fName = fieldNames{j};
            value = s.(fName);

            % Only write numeric fields
            if isnumeric(value)
                fprintf(fid, '--- %s.%s ---\n', sName, fName);
                if isscalar(value)
                    fprintf(fid, '%g\n', value);
                elseif isvector(value)
                    fprintf(fid, '%g ', value);
                    fprintf(fid, '\n');
                elseif ismatrix(value)
                    for row = 1:size(value, 1)
                        fprintf(fid, '%g ', value(row, :));
                        fprintf(fid, '\n');
                    end
                else
                    fprintf(fid, '[Unsupported array size >2D]\n');
                end
            else
                fprintf(fid, '--- %s.%s ---\n[Non-numeric type skipped]\n', sName, fName);
            end
        end
        fprintf(fid, '\n');  % Blank line between structs
    end

    fclose(fid);
    fprintf('Data written to ASCII file: %s\n', filename);
end