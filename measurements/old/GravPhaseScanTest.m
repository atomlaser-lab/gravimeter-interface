function GravPhaseScanTest(r)

if r.isInit()
    r.data.phase = 0:5:180;
    r.data.order = [-1,1];
    r.c.setup('var',r.data.order,r.data.phase);
elseif r.isSet()
    r.make(0,216.5e-3,1.05,0.2,r.data.phase(r.c(2)),5e-3,r.data.order(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, Order: %d, Phase: %.2f\n',r.c.now,r.c.total,r.data.order(r.c(1)),r.data.phase(r.c(2)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    nn = r.c.now;
    pause(0.1);
    c = Abs_Analysis_NClouds('last');
    if ~c(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    
    r.data.files{nn,1} = c(1).raw.files;

    r.data.N{i2,i1} = c.get('N');
    r.data.R(i2,i1) = c(2).N./sum(c.get('N'));
    r.data.Rsum(i2,i1) = c(2).Nsum./sum(c.get('Nsum'));
    
    if nn > 2
        figure(10);clf;
        for mm = 1:size(r.data.R,2)
            errorbar(r.data.phase(1:i2),r.data.R(1:i2,mm),0.005*ones(size(r.data.R(1:i2,mm))),'o-');
            hold on
            errorbar(r.data.phase(1:i2),r.data.Rsum(1:i2,mm),0.005*ones(size(r.data.Rsum(1:i2,mm))),'sq--');
        end
    end
    
end