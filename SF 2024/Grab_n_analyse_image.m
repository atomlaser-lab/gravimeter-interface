function Grab_n_analyse_image(imageNumber)
    % Check if an image number is provided, otherwise default to 1
    if nargin < 1
        imageNumber = 1; % Default to 1 if no input is provided
    end

    % Grab the image data for the specified image number
    img_details = grab_image_data_secondinstance('last', imageNumber); % Use the provided image number
    y_filt = img_details; % Use the appropriate data from img_details for fitting

    % Define fitting parameters
    spatial_freq = 0.1555; % Set the spatial frequency

    % Define the fitting range (you can modify based on the data size)
    beg = 10;
    finish = 2000;
    x_data = beg:finish;

    % Define the single Gaussian modulated by sine wave model
    single_gauss_mod_sine = @(b, x) ...
        b(1) * exp(-((x - b(2)).^2) / b(3)^2) .* (1 - b(4) * sin(spatial_freq * x - b(5))) + b(6); % Offset included

    % Initial parameters for the fit
    beta0_single = [13, 1500, 80, 0.8, -2, 5]; % Initial guess

    % Setting up the fitting options
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';
    opts.TolFun = 1e-19;
    opts.TolX = 1e-19;
    opts.Robust = 'on';
    opts.Display = 'off'; % Set to 'iter' for more details during the fitting
    opts.MaxIter = 1000;
    opts.UseParallel = true;

    % Perform the fitting
    try
        [beta_single, R, J, CovB, MSE] = nlinfit(x_data, y_filt(beg:finish), single_gauss_mod_sine, beta0_single, opts);
    catch ME
        disp(['Fitting failed: ' ME.message]);
        beta_single = NaN(1, length(beta0_single));
    end

    % Plotting the data and the fit
    figure(100); % Create a new figure or use figure 100
    clf; % Clear the figure
    hold on;
    plot(x_data, y_filt(beg:finish), 'b', 'LineWidth', 1.5); % Plot the filtered data
    plot(x_data, single_gauss_mod_sine(beta_single, x_data), 'r', 'LineWidth', 1.5); % Plot the fitted curve
    title(['Data and Fit for Image Number ' num2str(imageNumber)]);
    xlabel('Position');
    ylabel('Intensity');
    legend('Filtered Data', 'Fit');
    grid on;
    set(gca, 'FontSize', 14); % Improve readability
    hold off;

    % Print out the fitted parameters
    disp(['Fitted parameters for Image Number ' num2str(imageNumber) ':']);
    disp(beta_single);

    % Create a report figure to display the fitted parameters
    figure(101); % Create a new figure or use figure 101
    clf; % Clear the figure
    reportStr = sprintf([
        'Fitted Parameters for Image Number %d:\n\n', ...
        'Amplitude (b1): %.3f\n', ...
        'Center Position (b2): %.3f\n', ...
        'Width (b3): %.3f\n', ...
        'Contrast (Amplitude Modulation, b4): %.3f\n', ...
        'Phase (b5): %.3f\n', ...
        'Offset (b6): %.3f\n' ...
        ], imageNumber, beta_single(1), beta_single(2), beta_single(3), beta_single(4), beta_single(5), beta_single(6));
    
    % Display the report as text in the figure
    annotation('textbox', [0.1, 0.1, 0.8, 0.8], 'String', reportStr, 'FontSize', 14, 'EdgeColor', 'none');
    title(['Parameter Report for Image Number ' num2str(imageNumber)]);
    axis off;

    % Optionally, highlight key parameters, such as contrast, in the console
    disp(['Contrast (Amplitude Modulation, b4) for Image Number ' num2str(imageNumber) ': ' num2str(beta_single(4))]);
end