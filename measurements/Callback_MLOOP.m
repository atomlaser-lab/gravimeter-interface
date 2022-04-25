function Callback_MLOOP(r)

if r.isInit()
    %
    % Define some constants
    %
    r.data.constants.dir = 'C:\Users\admin\Desktop\matlab-control\gravimeter-interface';
    r.data.constants.input_file = fullfile(r.data.constants.dir,'exp_input.mat');
    r.data.constants.output_file = fullfile(r.data.constants.dir,'exp_output.mat');
    delete(r.data.constants.input_file);
    delete(r.data.constants.output_file);
    %
    % Set default parameters as a vector of values
    %
    r.data.default_params = 20e-3;
    %
    % Run forever
    %
    r.c.setup(Inf);
elseif r.isSet()
    %
    % Read input file.  Run only starts when there is an input file
    %
    max_read_attempts = 10;
    read_wait_time = 1;
    mm = 1;
    while ~isfile(r.data.constants.input_file)
        if mm > max_read_attempts
            fprintf('Experimental input file could not be found after %d attempts. Exiting',max_read_attempts);
            r.stop;
            return
        else
            pause(read_wait_time);
            mm = mm + 1;
        end
    end
    %
    % Read input file.  Assumed to be a .mat file
    %
    new_data = load(r.data.constants.input_file);
    new_params = new_data.params;
    % Convert to cell array
    full_param_list = [r.data.default_params(:);new_params(:)];
    args = mat2cell(full_param_list,ones(numel(full_param_list),1));
    %
    % Upload and print parameters
    %
    r.make(args{:}).upload;
    fprintf('Run %d/%d, Parameters =',r.c.now,r.c.total);
    fprintf(' %.5f',new_params);
    fprintf('\n');
    %
    % Store new parameters, delete input file
    %
    r.data.params(r.c.now(),:) = new_params;
    delete(r.data.constants.input_file);

elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.25);
    img = Abs_Analysis('last',1);
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and tells M-LOOP that the run was bad
        %
        cost = 10;
        uncer = 0;
        bad = true;
        save(r.data.constants.output_file,'cost','uncer','bad');
        return
    end
    %
    % Get data
    %
    r.data.files{i1,1} = img.raw.files;
    r.data.N(i1,1) = img.get('N');
    r.data.T(i1,1) = sqrt(prod(squeeze(img.get('T'))));
    r.data.PSD(i1,1) = img.get('PSD');
    r.data.OD(i1,1) = img.get('peakOD');
    %
    % Plot data
    %
    figure(124);clf;
    subplot(1,2,1);
    errorbar(1:i1,r.data.N(1:i1),0.05*r.data.N(1:i1),'o');
    plot_format('Run','Number','',12);
    ylim([0,Inf]);
    grid on;
    subplot(1,2,2);
    errorbar(1:i1,r.data.PSD(1:i1),0.05*r.data.PSD(1:i1),'o');
    plot_format('Run','PSD','',12);
    ylim([0,1e-4]);
    grid on;
    %
    % Write new output file. Remember - M-LOOP MINIMIZES the cost, so it
    % should be the negative of any positive measured quantity
    %
    if isnan(r.data.PSD(i1,1))
        cost = 10;
        uncer = 0;
        bad = false;
    elseif r.data.N(i1,1) < 1e6 || r.data.T(i1,1) > 30e-6
        cost = 10;
        uncer = 0;
        bad = false;
    else
        cost = -r.data.PSD(i1,1)*1e5;
        uncer = 0.05*r.data.PSD(i1,1)*1e5;
        bad = false;
    end
    save(r.data.constants.output_file,'cost','uncer','bad');
    r.data.cost(i1,:) = [cost,uncer];
    figure(83);clf;
    errorbar(1:i1,r.data.cost(:,1),r.data.cost(:,2),'o');
    grid on;
    ylim([-5,0]);
    plot_format('Run','Cost','',12);
end


end