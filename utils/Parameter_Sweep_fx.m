function Parameter_Sweep_fx(param2sweep,sweep_param_ar, lam0,crystal,I0,res,N,options, SaveDir)

Var_field = strsplit(param2sweep,'.'); % find variable name and field name to sweep from
StructName = Var_field(1);
FieldName = Var_field(2);

% Find exact match
Sweepable_variables.lam0 = lam0;
Sweepable_variables.crystal = crystal;
Sweepable_variables.I0 = I0;
Sweepable_variables.res = res;
filename = strcat('Sweep_param_',param2sweep);
% Create folder if it doesn't exist
if ~exist( SaveDir, 'dir')
mkdir( SaveDir);
end
% Write all variables to ASCII 
write_allStructs_to_ascii(Sweepable_variables, fullfile(SaveDir,filename))

for i = 1:length(sweep_param_ar)
    if isfield(Sweepable_variables, StructName)
        sweep_var_struct = Sweepable_variables.(StructName{1});
        if isfield(sweep_var_struct,FieldName)
            
            Sweepable_variables.(StructName{1}).(FieldName{1}) = sweep_param_ar(i);
        else
            error('The struct field is not in the list of sweepable variables See Parameter_Sweep_fx.m')
   
        end


    else
        error('The struct is not in the list of sweepable variables See Parameter_Sweep_fx.m')
    end
    % Run parameter i
    filename_i = strcat(filename,'_',num2str(sweep_param_ar(i)));
    ff = figure(2); % Set figure to two so we don't overload GPU with grpahics
    fig_num = ff.Number;
    TROPO_Generalized(Sweepable_variables.lam0,...
        Sweepable_variables.crystal,...
        Sweepable_variables.I0,...
        Sweepable_variables.res,...
        N)
    ff = gcf;
    run('tilefigs.m')
    % Take Screenshot of all the figures and then save to director
    % Takes a screenshot and saves it to the specified folder


    % Take screenshot of entire screen
    drawnow
    pause(1); % Pause so that it can load the graphics before screenshot of screen
    img = screencapture(0);

    % Convert to image (if needed)
    % screenshot = frame2im(img);
    imwrite(img, fullfile(SaveDir,strcat(filename_i,'.png'))); %save too dir


end
    

end



