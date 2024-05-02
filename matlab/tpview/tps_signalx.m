classdef tps_signalx < hgsetget
    properties
        tag
        active = true;
        sel         
        spikepar
        kcond            % condition number; 0=not a simple condition 
        delay = 0;       % depends only on selection
        delayshown = 0;  % depends on opdef and selection (=delay if show delay, =0 otherwise)
        tidx = [];       % depends on opdef and selection (if show delay) 
        data = [];       % depends on data, datamode and selection
        data2 = [];
        dataop = [];     % depends on everything, even on other selections!
        data2op = [];     
        spikes = [];
        spikefit = [];
        validspike = false % in many situations the data on which spikes are calculated is changed, but spikes are not deleted, but they are not 'valid' 
        spikes2 = [];
        spikefit2 = [];
        validspike2 = false % in many situations the data on which spikes are calculated is changed, but spikes are not deleted, but they are not 'valid' 
    end
    methods
        function x2 = copy(x)
            if ~isscalar(x)
                x2 = x;
                for k=1:numel(x), x2(k)=copy(x(k)); end
                return
            end
            x2 = tps_signalx;
            x2.sel = x.sel;
            x2.tidx = x.tidx;
            x2.data = x.data;
            x2.data2 = x.data2;
            x2.dataop = x.dataop;
            x2.data2op = x.data2op;
            x2.spikes = x.spikes;
            x2.spikefit = x.spikefit;
            x2.validspike = x.validspike;
            x2.spikes2 = x.spikes2;
            x2.spikefit2 = x.spikefit2;
            x2.validspike2 = x.validspike2;
        end
    end
end