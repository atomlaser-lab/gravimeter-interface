classdef DDSChannel < TimingControllerChannel
    %DDSChannel Defines a DDS channel as an extension
    %to the TimingControllerChannel class.  DDSChannel objects have values
    %that have 3 sub-values for frequency, normalized power, and phase
    
    properties
        channel
        calibrationData
        rfscale
        power_conversion_method
    end
    
    properties(Constant)
        CLK = 1e6;
        DEFAULT_FREQ = 110;
        POWER_CONVERSION_MODEL = 'model';
        POWER_CONVERSION_DBM_INTERP = 'dbm_interp';
        POWER_CONVERSION_HEX_INTERP = 'hex_interp';
    end
    
    methods
        function ch = DDSChannel
            ch = ch@TimingControllerChannel;
            ch.default = [110,0,0];
            ch.bounds = [80,    0,  0;
                         150,   1,  1200];
            ch.IS_DDS = true;
            ch.power_conversion_method = 'model';
        end
        
        function ch = setDefault(ch,v)
            if numel(v) == 1 && v == 0
                v = [110,0,0];
            end
            ch = setDefault@TimingControllerChannel(ch,v);
        end
        
        function ch = checkValue(ch,v)
            %CHECKVALUE Checks a given value, given as freq, amp, and
            %phase, is within the given ranges
            %
            %   CH = CHECKVALUE(CH,V) checks a given value V, given as
            %   freq, power, and phase, is within the given ranges
            str = {'Frequency','Power','Phase'};
            for nn = 1:size(v,2)
                if any(v(:,nn) < ch.bounds(1,nn)) || any(v(:,nn) > ch.bounds(2,nn))
                    idx = find((v(:,nn) < ch.bounds(1,nn)) | (v(:,nn) > ch.bounds(2,nn)),1,'first');
                    error('%s %f is out of range [%f,%f]',str{nn},v(idx,nn),ch.bounds(1,nn),ch.bounds(2,nn));
                end
            end
        end
        
        function ch = expand(ch,tin)
            ch.sort;
            tin = tin(:);
            fnew = interp1(ch.times,ch.values(:,1),tin,'previous','extrap');
            pownew = interp1(ch.times,ch.values(:,2),tin,'previous','extrap');
            phnew = interp1(ch.times,ch.values(:,3),tin,'previous','extrap');
            ch.times = tin;
            ch.values = [fnew,pownew,phnew];
            ch.numValues = numel(ch.times);
        end
        
        function data = compile(ch,delay)
            data.power_conversion_method = ch.power_conversion_method;
            if ch.numValues == 0
                data.t = 0;
                data.freq = ch.default(1);
                data.pow = ch.default(2);
                data.phase = ch.default(3);
                return
            end
            t = ch.times - delay;
            idx = (t >= 0);
            t = t(idx);
            v = ch.values(idx,:);
            
            if isempty(t)
                t = 0;
                v = [ch.DEFAULT_FREQ,0,0];
            elseif t(1) ~= 0
                t = [0;t];
                v = [ch.DEFAULT_FREQ,0,0;v];
            end
            data.t = t;
            data.freq = v(:,1);
            if strcmpi(ch.power_conversion_method,ch.POWER_CONVERSION_MODEL)
                data.pow = 30 + 10*log10(ch.opticalToRF(v(:,2),1,ch.rfscale));
            elseif strcmpi(ch.power_conversion_method,ch.POWER_CONVERSION_DBM_INTERP)
                data.pow = 30 + 10*log10(ch.opticalToRF(v(:,2),ch.calibrationData));
            elseif strcmpi(ch.power_conversion_method,ch.POWER_CONVERSION_HEX_INTERP)
                data.pow = ch.opticalToHex(v(:,2),ch.calibrationData);
            else
                error('Need to supply only one of calibration data or RF scale to convert optical power to RF power!');
            end
            data.phase = v(:,3);
        end
        
        function [tplot,vplot] = getPlotValues(ch,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                offset = 0;
                finalTime = [];
                returnHandle = false;
                plotIdx = 1:size(ch.values,2);
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case 'offset'
                            offset = v;
                        case 'finaltime'
                            finalTime = v;
                        case 'returnhandle'
                            returnHandle = v;
                        case 'plotidx'
                            plotIdx = v;
                    end
                end
            end
            [t,v] = ch.getEvents;
            if ~ch.exists && ~returnHandle
                tplot = [];
                vplot = [];
                return
            end
            
            if t(end) ~= finalTime
                t = [t;finalTime];
                v = [v;v(end,:)];
            end
%             tplot = sort([t;t-1/TimingSequence.SAMPLE_CLK]);
%             tplot = tplot(tplot >= 0);
%             vplot = zeros(numel(tplot),numel(plotIdx));
%             for nn = 1:numel(plotIdx)
%                 vplot(:,nn) = interp1(t,v(:,plotIdx(nn)),tplot,'previous');
%             end
            tplot = t;
            vplot = v(:,plotIdx);
            vplot = vplot + offset;
        end

    end
    
    methods(Static)
        function rf = opticalToRF(P,varargin)
            if numel(varargin) == 1 && isstruct(varargin{1})
                data = varargin{1};

                if isfield(data,'Pmax') == 0
                    data.Pmax = max(data.optical_power);
                end
                P = P*data.Pmax;
                if isfield(data,"Popt")
                    rf = interp1(data.Popt,data.Prf,P,'pchip');
                else
                    if isfield(data,'amp') == 0 && isfield(data,'amp_dBm') == 1
                        data.amp = (10.^(data.amp_dBm/10))*10^-3;
                    end
                    rf = interp1(data.optical_power,data.amp,P,'pchip');
                end
                % manually set zeros
                rf(P==0) = 1e-6;

            elseif numel(varargin) == 2
                Pmax = varargin{1};
                rfmax = varargin{2};
                rf = (asin((P/Pmax).^0.25)*2/pi).^2*rfmax;
            end
        end
        
        function amp = opticalToHex(P,varargin)
            data = varargin{1};
            [~,k] = sort(data.amp);
            data.amp = data.amp(k);data.optical_power = data.optical_power(k);
            [max_optical_power,idx] = max(data.optical_power);
            data.amp = data.amp(1:idx);data.optical_power = data.optical_power(1:idx);
            P = P.*max_optical_power;
            amp = interp1(data.optical_power,data.amp,P,'pchip');
        end
    end
    
end