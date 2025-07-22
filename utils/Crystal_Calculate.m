

function [crystal,lam0] = Crystal_Calculate(crystal,lam0)
% Define some values to use in specifying set of inputs for calling snlo_qmix_func.
% These values will be stuck in the data structure below.
QPM = crystal.QPM; % yes or no for QPM 
PM_type = crystal.PM_type; % type-0, type-1, or type-2
crystal_name = crystal.name;              % name of crystal. string. must match an entry in the file crystal_list.txt.
temperature = crystal.temperature;                  % temperature in kelvin for use in temperature dependent Sellmeier calculations.
% global c eps_0;  
c =  2.99792458E8;           % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;           % dielectric constant in vacuum, F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(crystal.QPM,'no')
    plane = crystal.plane;                       % plane (ignored for uniaxial crystals); should be a pair of characters 'XY', 'XZ', or 'YZ'
    lambdap = lam0.p*1e9; % Must give a pump wavlength
    if ~isfield(lam0,'i') % no idler defined
        lambdas = lam0.s*1e9;
        lambdai = 0;
    else
        lambdai = lam0.i*1e9;
        lambdas = 0;
    end
    wavelengths = [lambdai,lambdas,lambdap];          % wavelengths (in nm) (here one must be zero)
    type = 'Mix';                       % mixing type ('Mix' or 'OPO')

    
    outputfilename = 'Qmix.dat';
    if strcmpi(crystal.PM_type,'type-0')
        warning('Crystal Phase matching choice of type-0 is not valid for non-QPM crystal')
    end
    % use Qmix to calculate crystal parameters
    % construct data structure to use as input when calling snlo_qmix_func
    % input_set.qmix_selected_crystal = crystal_name;
    input_set.qmix_temperature = temperature;
    input_set.qmix_principal_plane = plane;
    input_set.qmix_wavelength_red1 = wavelengths(1); 
    input_set.qmix_wavelength_red2 = wavelengths(2); 
    input_set.qmix_wavelength_blue = wavelengths(3); 
    input_set.qmix_type = type;
    input_set.qmix_crystal_name = crystal_name;


    crystal_list = importdata('crystal_list.txt'); % load list of crystals
    [~,input_set.qmix_selected_crystal] = max(strcmpi(crystal_list,input_set.qmix_crystal_name));
    problem = snlo_qmix_func(input_set);
    run_function = problem{1};
    % simulate pressing the 'Run' button
    run_function();


    fid = fopen(outputfilename, 'r');
    filestr = fscanf(fid,'%s'); % load entire contents into variable filestr 
    num_phasematches = length(strfind(filestr,'Walkoff')); % find length of matches for the word 'Walkoff' in contents of Qmix.dat
    frewind(fid); % go back to the beginning of 


    K = 1;
%     % stick some appropriately sized arrays of 0 in the output cell arrays
%     output_wavelengths{K} = zeros(3,num_phasematches);          % 3 entries per phase match: wavelength of red1, red2, blue
%     output_polarizations{K} = zeros(3,num_phasematches);        % 3 entries: polarizations (o or e) for red1, red2, and blue
%     output_walkoffs{K} = zeros(3,num_phasematches);             % 3 entries: walk off of red1, red2, blue
%     output_phase_velocities{K} = zeros(3,num_phasematches);     % 3 entries: phase velocity indices of red1, red2, blue
%     output_group_velocities{K} = zeros(3,num_phasematches);     % 3 entries: group velocity indices of red1, red2, blue
%     output_gdd{K} = zeros(3,num_phasematches);                  % 3 entries: group delay dispersions of red1, red2, blue
%     output_theta{K} = zeros(1,num_phasematches);                % 1 entry: phase match angle
%     output_d_eff{K} = zeros(1,num_phasematches);                % 1 entry: d effective
%     output_SoL2{K} = zeros(1,num_phasematches);                 % 1 entry: S_0 * L^2
%     output_angle_tol{K} = zeros(1,num_phasematches);            % angle tolerance
%     output_temp_range{K} = zeros(1,num_phasematches);           % temperature range
%     output_mix_accpt_ang{K} = zeros(2,num_phasematches);        % acceptance angles
%     output_mix_accpt_bw{K} = zeros(2,num_phasematches); 

    % loop through each phase match for this temperature
    for K = 1:num_phasematches

        line_contents = fgetl(fid);
        str_new = regexprep(line_contents,'[^a-zA-Z\s]','');
        output_PM_type{K} = str_new(str_new~=' ');
        % first line of Qmix output: for each wave, wavelength and polarization listed in
        parsed_contents = sscanf(line_contents, '%f(%c) +  %f(%c)  =  %f(%c)');
        output_wavelengths{K} = [parsed_contents(1), parsed_contents(3), parsed_contents(5)];
        output_polarizations{K} = char([(parsed_contents(2)),(parsed_contents(4)),...
            (parsed_contents(6))]);
        
        % second line: walk off angles
        line_contents = fgetl(fid);
        output_walkoffs{K} = sscanf(line_contents, 'Walkoff [mrad]      = %f %f %f');

        % third line: phase velocity indices
        line_contents = fgetl(fid);
        output_phase_velocities{K} = sscanf(line_contents, 'Phase velocities    = c/ %f %f %f');

        % third line: group velocity indices
        line_contents = fgetl(fid);
        output_group_velocities{K} = sscanf(line_contents, 'Group velocities    = c/ %f %f %f');

        % fourth line: group delay dispersions
        line_contents = fgetl(fid);
        output_gdd{K} = sscanf(line_contents, 'GrpDelDisp(fs^2/mm) = %f %f %f');

        % fifth line: phase match angle
        line_contents = fgetl(fid);
        output_theta{K} = sscanf(line_contents, 'At theta            = %f deg.');

        % sixth line: d effective
        line_contents = fgetl(fid);
        output_d_eff{K} = sscanf(line_contents, 'd_eff               = %f pm/V');

        % seventh line: product of characteristic irradiance S_o and square of crystal length
        line_contents = fgetl(fid);
        output_SoL2{K} = sscanf(line_contents, 'S_o * L^2           = %f watt');

        % eigth line: crystal angle tolerance
        line_contents = fgetl(fid);
        output_angle_tol{K} = sscanf(line_contents, 'Crystal ang. tol.   = %f mrad-cm');

        % ninth line: temperature range
        line_contents = fgetl(fid);
        output_temp_range{K} = sscanf(line_contents, 'Temperature range   = %f K-cm');

        % mix accptance angles *** note: for OPO mix type rather than 'Mix' mix type, first 3 characters of sscanf are OPO rather than Mix
        line_contents = fgetl(fid);
        output_mix_accpt_ang{K} = sscanf(line_contents, 'Mix accpt ang   = %f %f mrad-cm');

        % mix acceptance bandwidths *** note: for OPO mix type rather than 'Mix' mix type, first 3 characters of sscanf are OPO rather than Mix
        line_contents = fgetl(fid);
        output_mix_accpt_bw{K} = sscanf(line_contents, 'Mix accpt bw    = %f %f cm^-1-cm');

        % the last line will be blank
        line_contents = fgetl(fid);
    end
    fclose(fid); % close file handle
elseif strcmpi(crystal.QPM,'yes')
    outputfilename = 'qpm.dat';
    lambdap = lam0.p*1e9; % Must give a pump wavlength
    if ~isfield(lam0,'i') % no idler defined
        lambdas = lam0.s*1e9;
        lambdai =1/(1/lambdap-1/lambdas);
        lam0.i = lambdai./1e9;
    else
        lambdai = lam0.i*1e9;
        lambdas = 1/(1/lambdap-1/lambdai);
        lam0.s = lambdas./1e9;
    end
    wavelengths = ([lambdai+10,lambdas-10,lambdap]);          % wavelengths (in nm) (here one must be zero)
    switch crystal.PM_type
        case 'type-0'
            lam0s_pol = 'Z-pol.';
            lam0i_pol = 'Z-pol.';
            lam0p_pol = 'Z-pol.';
            pol = 'ZZZ'
        case 'type-1'
            lam0s_pol = 'X-pol.';
            lam0i_pol = 'X-pol.';
            lam0p_pol = 'Z-pol.';
            pol = 'XXZ'
        case 'type-2'
            lam0s_pol = 'X-pol.';
            lam0i_pol = 'Z-pol.';
            lam0p_pol = 'Z-pol.';
            pol = 'XZZ'
        
    end
    qpm_crystal_list = importdata('crystal_list_qpm.txt'); % load list of crystals
    inputs.qpm_temperature = crystal.temperature;               % temperature (in K)
   
    inputs.qpm_wavelengths = sort(wavelengths);    % pump, bluest, reddest wavelengths [nm]
    inputs.currentCrystalPopupName = crystal.name;   % specify one of the strings which are contained in qpm crystal list (crystal_list_qpm.txt)
    inputs.qpm_currentCrystalPopupValue = find_selection_number_from_string(qpm_crystal_list,inputs.currentCrystalPopupName); % find the position in the list of crystals which corresponds to the name of the crystal you specified
        
    inputs.qpm_currentPolValue = crystal.qpm_currentPolValue;             % position in list of polarizations which corresponds to the one you want to use (run qpm and select your crystal to find the list of polarizations, and which set this value to the position in the list of polarizations)
        a = snlo_qpm_func(inputs);
    inputs.qpm_temptune_temp_range = [273,450]; % temperature range for temperature tuning [K] - not used in this example
    inputs.qpm_temptune_period = 10;            % period for temperature tuning [um] - not used in this example
    inputs.qpm_pumptune_wavelength_range = [796,810]; % wavelength range for pump tuning [nm] - not used in this example
    inputs.qpm_pumptune_period = 28;            % period for pump tuning [um] - not used in this example

    
    % call qpm, run model, load output
    clear problem;
    outputs = cell(1); % pre-allocate a cell array to stick contents of model output file into
    K = 1;
        problem(K) = inputs;  % copy the inputs specified earlier as the starting point for this problem's inputs
        fcn_handles = snlo_qpm_func(problem(K)); % call the SNLO qpm function with the problem set, and assign the returned cell array of function handles which are local to that file to 'fcn_handles'
        % run poling calculation
        run_fcn = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
        outputs{K} = importdata('qpm.dat'); % load the model output from the file qpm.dat, stick it in an element in the cell array 'output'
            Poling = abs([outputs{K}(:,3);outputs{K}(end:-1:1,3)]);
            lam_pol = [outputs{K}(:,1);outputs{K}(end:-1:1,2)];
        [~,inx] = min(abs(lam_pol -lam0.s));
        crystal.Period = Poling(inx).*1e-6;
        f = gcf;
        Tit = f.Children.Title.String;
        deff =sscanf(Tit{2},['T = 300.00 K, Pol.(s,i,p)=',pol,', d_{eff} = %f pm/V']);
        % pol =sscanf(Tit{2},'T = 300.00 K, Pol.(s,i,p)=%c, d_{eff} = deff pm/V');
         
        % run gvm calculation
        run_fcn = fcn_handles{12}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
       
        axObj = f.Children;
        dataObj = axObj(strcmpi(get(axObj,'Type'),'axes')).Children;
        ng_pump = regexp(dataObj(strcmpi(get(dataObj,'Type'),'Text')).String,'\d*','Match');
        ng_pump = str2num(strcat(ng_pump{1},'.',ng_pump{2}));
        
        ng_signal = findobj(dataObj, 'Type', 'line','Displayname',lam0s_pol);
        if isempty(ng_signal)
            ng_signal = findobj(dataObj, 'Type', 'line');
            ng_signal = ng_signal(1);
        end
            [~,inx] = min(abs(ng_signal.XData - lambdas));
            ng_signal = ng_signal.YData(inx);
        ng_idler = findobj(dataObj, 'Type', 'line','Displayname',lam0s_pol);
        if isempty(ng_idler)
            ng_idler = findobj(dataObj, 'Type', 'line');
            ng_idler = ng_idler(1);
        end
            [~,inx] = min(abs(ng_idler.XData - lambdai));    
            ng_idler = ng_idler.YData(inx);
            
        run_fcn = fcn_handles{13}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
        run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
       
        axObj = f.Children;
        dataObj = axObj(strcmpi(get(axObj,'Type'),'axes')).Children;
        GDD_pump = regexp(dataObj(strcmpi(get(dataObj,'Type'),'Text')).String,'\d*','Match');
        GDD_pump = str2num(strcat(GDD_pump{1}));
        
            
        GDD_signal = findobj(dataObj, 'Type', 'line','Displayname',lam0s_pol);
        if isempty(GDD_signal)
            GDD_signal = findobj(dataObj, 'Type', 'line');
                mm1 = min(abs(GDD_signal(1).XData - lambdas));
                mm2 = min(abs(GDD_signal(2).XData - lambdas));
                [~,inx] =min([mm1,mm2]);
            GDD_signal = GDD_signal(inx);
        end    
            [~,inx] = min(abs(GDD_signal.XData - lambdas));
            GDD_signal = GDD_signal.YData(inx);
        GDD_idler = findobj(dataObj, 'Type', 'line','Displayname',lam0s_pol);
        if isempty(GDD_idler)
            GDD_idler = findobj(dataObj, 'Type', 'line');
                mm1 = min(abs(GDD_idler(1).XData - lambdai));
                mm2 = min(abs(GDD_idler(2).XData - lambdai));
                [~,inx] =min([mm1,mm2]);
            GDD_idler = GDD_idler(inx);
        end       
            [~,inx] = min(abs(GDD_idler.XData - lambdai));    
            GDD_idler = GDD_idler.YData(inx);
       
            
    % Save values
    crystal.beta1_p = (c./ng_pump).^-1; crystal.beta1_s = (c./ng_signal).^-1; crystal.beta1_i = (c./ng_idler).^-1;
    crystal.beta2_p = GDD_pump.*1e-15*1e-15./1e-3; crystal.beta2_s = GDD_signal.*1e-15*1e-15./1e-3; crystal.beta2_i = GDD_idler.*1e-15*1e-15./1e-3;
    crystal.deff = deff*1e-12;
    run_fcn = fcn_handles{1}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
    run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
    
    f = gcf;
        Tit = f.Children.Title.String;
    run_fcn = fcn_handles{13}; % the function handle of the 'Run' button callback is the first returned; the 'run' button callback function is always the first element of the cell array of returned function handles
    run_fcn(); % call the 'Run' button callback - equivalent to clicking the run button
    
    hold on;
    plot(lam0.i*1e9,GDD_idler,'*r')
    plot(lam0.s*1e9,GDD_signal,'*g')
    grid on;
    title({[strcat(num2str(round(lam0.i.*1e9)),'nm Idler','    |',num2str(round(lam0.s.*1e9)),'nm Signal  |',Tit{1})],Tit{2}})
    drawnow
end
    
    
end
