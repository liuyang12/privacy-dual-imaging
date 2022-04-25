function [meas] = readmeasdata(rawdata_path,params)
%READMEASDATA Read measurement data from the raw data.
%   See also READRAWPCIDATA.
includelastframe = 0;
if isfield(params,'includelastframe'), includelastframe = params.includelastframe; end
% [1] read raw PCI data
vall = readrawPCIdata(rawdata_path,params);

% [2] seperate the measure data for each channel 
for ich = 1:length(params.chnum)
    v = vall(:,params.chnum(ich));
    if max(abs(v))>max(v)
        v = -v;
    end
    % [2.1] extract the measurement between adjacent pattern switching
    % figure; plot(v,'b'); grid on; hold on; plot(params.measthresh*ones(size(v)),'r');
    ipat = find(v>params.measthresh); % indexes above measurement threshold
    dpat = find(diff(ipat)>2); % dramatic drop of adjacent switching
    % patno = [ipat([1;dpat+1])+50,ipat([dpat;end])-30]; % pattern indexes (from the first column to the second)
    patno = [ipat([1;dpat+1]),ipat([dpat;end])]; % pattern indexes (from the first column to the second)
    if ~includelastframe % exclude the last frame (for WinTech DMD)
        patno = patno(1:end-1,:); % exclude the last frame (excessive, 1st)
    end
    % meas(:,ich) = v(floor(0.3*patno(:,1)+0.7*patno(:,2))); % measurement for each channel
    t = 1;
    interval = round((patno(end,1)-patno(1,1)+patno(end,2)-patno(1,2))/length(patno)/2);
    patno(patno(:,2)-patno(:,1)<0.5*interval,:) = []; % exclude extra wrong data resulted from threshold
    interval = round((patno(end,1)-patno(1,1)+patno(end,2)-patno(1,2))/length(patno)/2);
    for p = 1:length(patno)
        if p>1
            p1 = patno(p,1);
            p0 = patno(p-1,1); q0 = patno(p-1,2);
            while p1-p0>1.5*interval
                p0 = p0+interval; q0 = q0+interval;
                meas(t,ich) = mean(v(floor(2/5*p0+3/5*q0):ceil(1/5*p0+4/5*q0)));
                t = t+1;
            end
        end     
        meas(t,ich) = mean(v(floor(2/5*patno(p,1)+3/5*patno(p,2)):ceil(1/5*patno(p,1)+4/5*patno(p,2))));
        t = t+1;
    end
end

end

