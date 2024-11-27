classdef FeedbackOptions < SequenceOptionsAbstract
    
    properties
        enable_ndi
        pulse_power
        pulse_time
        pulse_delay
        cycle_time

        num_images
        ref_images

        enable_fb_laser
        fb_laser_power
    end

    methods
        function self = FeedbackOptions(varargin)
            self.setDefaults;
            self = self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.enable_ndi = 1;
            self.pulse_power = 0.1;
            self.pulse_time = 5e-6;
            self.pulse_delay = 50e-6;
            self.cycle_time = 2e-3;
            self.num_images = 0;
            self.ref_images = 2;
            self.enable_fb_laser = 0;
            self.fb_laser_power = 0;
        end

    end
end