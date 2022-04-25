function [meas] = readvidmeasdata(I,params)
%READMEASDATA Read measurement data from the raw data.
%   See also READRAWPCIDATA.
BLANKSEG = true; % measurement segmentation using blank frames (sharp drops)
includelastframe = 1; % flag of including last frame or not
avginterval = [0.3, 0.8]; % the interval to be averaged for each pattern
interfactor = 1.75; % interval factor (to determine the segmentation)
if isfield(params,'BLANKSEG'),                 BLANKSEG = params.BLANKSEG;      end
if isfield(params,'includelastframe'), includelastframe = params.includelastframe; end
if isfield(params,'avginterval'),           avginterval = params.avginterval;      end
if isfield(params,'interfactor'),           interfactor = params.interfactor;      end

left  = avginterval(1);
right = avginterval(2);

% % [1] read raw PCI data
% vall = readrawPCIdata(rawdata_path,params);
[L, nch] = size(I); % length and number of color channels of the signal
I_total = sum(I,2); % seperation based on the total intensity
% [2] seperate the measure data for each channel 
% for ich = 1:length(params.chnum)
    % I = vall(:,params.chnum(ich));
    % [2.1] extract the measurement between adjacent pattern switching
    figure; plot(I_total,'b'); grid on; hold on; plot(params.measthresh*ones(size(I)),'r'); ylabel('lux');
    ipat = find(I_total>params.measthresh); % indexes above measurement threshold
    
    if ~BLANKSEG % segmentation without blank frames based on timing and total number of frames
        M = params.nmeas;
        t0 = ipat(1);
        tM = ipat(end);
        indall = round(t0+((1:M)-1/2)*(tM-t0)/M);
        % meas = I(indall,:)';
        meas = ((I(indall,:)+I(indall+1,:))/2)';
    else % segmentation with blank frames based on sharp drops of blank frames
        dpat = find(diff(ipat)>2); % dramatic drop of adjacent switching
        % patno = [ipat([1;dpat+1])+50,ipat([dpat;end])-30]; % pattern indexes (from the first column to the second)
        patno = [ipat([1;dpat+1]),ipat([dpat;end])]; % pattern indexes (from the first column to the second)
        if ~includelastframe % exclude the last frame (for WinTech DMD)
            patno = patno(1:end-1,:); % exclude the last frame (excessive, 1st)
        end
        % meas(:,ich) = v(floor(0.3*patno(:,1)+0.7*patno(:,2))); % measurement for each channel
        t = 1;
        interval = round((patno(end,1)-patno(1,1)+patno(end,2)-patno(1,2))/length(patno)/2);
        patno(patno(:,2)-patno(:,1)<0.2*interval,:) = []; % exclude extra wrong data resulted from threshold
        interval = round((patno(end,1)-patno(1,1)+patno(end,2)-patno(1,2))/length(patno)/2);
        for p = 1:length(patno)
            if p>1
                p1 = patno(p,1);
                p0 = patno(p-1,1); q0 = patno(p-1,2);
                while p1-p0>interfactor*interval
                    p0 = p0+interval; q0 = q0+interval;
                    for ich = 1:nch
                        meas(ich,t) = mean(I(floor((1-left)*p0+left*q0):ceil((1-right)*p0+right*q0),ich));
                    end
                    t = t+1;
                end
            end     
            for ich = 1:nch
                meas(ich,t) = mean(I(floor((1-left)*patno(p,1)+left*patno(p,2)):ceil((1-right)*patno(p,1)+right*patno(p,2)),ich));
            end
            t = t+1;
        end
    end
end

% end

