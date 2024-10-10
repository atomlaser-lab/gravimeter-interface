function Grab_n_analyse_dualSF(imageNumber)
    % Set the directory path where the image data is located
    directory = 'Z:\Gravy\2024\last\Initial velocity\Initial velocity combo interfereomter 1ms [ perhaps retake]';

    % Check if an image number is provided, otherwise default to 1
    if nargin < 1
        imageNumber = 1; % Default to 1 if no input is provided
    end

    % Grab the image data for the specified image number from the directory
    img_details = grab_image_data_secondinstance(directory, 'last', imageNumber);
    y_filt = img_details; % Use the appropriate data from img_details for fitting

    % Define fitting parameters
    spatial_freq = 0.1555; % Set the spatial frequency

    % Define the fitting range (you can modify based on the data size)
    beg = 10;
    finish = 2000;
    x_data = beg:finish;

    % Define the dual Gaussian modulated by sine wave model
    dual_gauss_mod_sine = @(b, x) ...
        (b(1) * exp(-((x - b(2)).^2) / b(3)^2) .* (1 - b(4) * sin(spatial_freq * x - b(5))) + ...
         b(6) * exp(-((x - b(7)).^2) / b(8)^2) .* (1 - b(9) * sin(spatial_freq * x - b(10))) + ...
         b(11)); % Offset included

    % Initial parameters for the fit
    beta0_dual = [13, 680, 80, 0.8, -2, 21, 1293, 80, 0.6, 0, 5]; % Initial guess for dual Gaussian

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
        [beta_dual, R, J, CovB, MSE] = nlinfit(x_data, y_filt(beg:finish), dual_gauss_mod_sine, beta0_dual, opts);
    catch ME
        disp(['Fitting failed: ' ME.message]);
        beta_dual = NaN(1, length(beta0_dual));
    end

    % Plotting the data and the fit
    figure(100); % Create a new figure or use figure 100
    clf; % Clear the figure
    hold on;
    plot(x_data, y_filt(beg:finish), 'b', 'LineWidth', 1.5); % Plot the filtered data
    plot(x_data, dual_gauss_mod_sine(beta_dual, x_data), 'r', 'LineWidth', 1.5); % Plot the fitted curve
    title(['Data and Fit for Image Number ' num2str(imageNumber)]);
    xlabel('Position');
    ylabel('Intensity');
    legend('Filtered Data', 'Fit');
    grid on;
    set(gca, 'FontSize', 14); % Improve readability
    hold off;

    % Print out the fitted parameters
    disp(['Fitted parameters for Image Number ' num2str(imageNumber) ':']);
    disp(beta_dual);

    % Create a report figure to display the fitted parameters
    figure(101); % Create a new figure or use figure 101
    clf; % Clear the figure
    reportStr = sprintf([
        'Fitted Parameters for Image Number %d:\n\n', ...
        'Gaussian 1:\n', ...
        '  Amplitude (b1): %.3f\n', ...
        '  Center Position (b2): %.3f\n', ...
        '  Width (b3): %.3f\n', ...
        '  Contrast (Amplitude Modulation, b4): %.3f\n', ...
        '  Phase (b5): %.3f\n\n', ...
        'Gaussian 2:\n', ...
        '  Amplitude (b6): %.3f\n', ...
        '  Center Position (b7): %.3f\n', ...
        '  Width (b8): %.3f\n', ...
        '  Contrast (Amplitude Modulation, b9): %.3f\n', ...
        '  Phase (b10): %.3f\n\n', ...
        'Offset (b11): %.3f\n' ...
        ], imageNumber, beta_dual(1), beta_dual(2), beta_dual(3), beta_dual(4), beta_dual(5), ...
        beta_dual(6), beta_dual(7), beta_dual(8), beta_dual(9), beta_dual(10), beta_dual(11));
    
    % Display the report as text in the figure
    annotation('textbox', [0.1, 0.1, 0.8, 0.8], 'String', reportStr, 'FontSize', 14, 'EdgeColor', 'none');
    title(['Parameter Report for Image Number ' num2str(imageNumber)]);
    axis off;

    % Optionally, highlight key parameters, such as contrast, in the console
    disp(['Contrast (Amplitude Modulation, b4 & b9) for Image Number ' num2str(imageNumber) ':']);
    disp(['  Contrast 1: ' num2str(beta_dual(4))]);
    disp(['  Contrast 2: ' num2str(beta_dual(9))]);
end
