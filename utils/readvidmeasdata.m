function [measall, patsegall] = readvidmeasdata(I,params)
%READMEASDATA Read measurement data from the raw data.
%   See also READRAWPCIDATA.
measthresh = 'median'; % median value of the raw measurement as the segement threshold
BLANKSEG = true; % measurement segmentation using blank frames (sharp drops)
includelastframe = 1; % flag of including last frame or not
avginterval = [0.3, 0.8]; % the interval to be averaged for each pattern
interfactor = 1.75; % interval factor (to determine the segmentation per pattern measurement)
segimfactor = 3; % per-image segmentation factor of the data point interval
processmeth = 'average'; % processing method to obtain each measurement
if isfield(params,'measthresh'),             measthresh = params.measthresh;       end
if isfield(params,'BLANKSEG'),                 BLANKSEG = params.BLANKSEG;         end
if isfield(params,'includelastframe'), includelastframe = params.includelastframe; end
if isfield(params,'avginterval'),           avginterval = params.avginterval;      end
if isfield(params,'interfactor'),           interfactor = params.interfactor;      end
if isfield(params,'processmeth'),           processmeth = params.processmeth;      end

left  = avginterval(1);
right = avginterval(2);

% % [1] read raw PCI data
% vall = readrawPCIdata(rawdata_path,params);
repeat = 1;
I_tmp = repmat(I',[repeat,1]);
I = I_tmp(:);
[L, nch] = size(I); % length and number of color channels of the signal
I_total = sum(I,2); % seperation based on the total intensity
if ischar(measthresh)
    switch lower(measthresh)
        case 'median'
            measthresh = median(I_total);
        case 'mean'
            measthresh = mean(I_total);
        case {'mixed', 'mix'}
            measthresh = (median(I_total) + mean(I_total)) / 2;
        otherwise
            disp('Unsupported thresholding method. Using MEDIAN instead.')
            measthresh = median(I_total);
    end
end
fprintf('  Segement threshold %.1f lux out of max %d lux. \n', measthresh, max(I_total)); 

% [2] seperate the measure data for each channel 
% for ich = 1:length(params.chnum)
    % I = vall(:,params.chnum(ich));
    % [2.1] extract the measurement between adjacent pattern switching
    % figure; plot(I_total,'b'); grid on; hold on; plot(measthresh*ones(size(I_total)),'r'); ylabel('lux');
    % ipat = find(I_total>measthresh); irange=ipat(1):ipat(end);figure; plot((1:length(irange))/length(irange),I_total(irange),'b'); grid on; hold on; plot((1:length(irange))/length(irange),measthresh*ones(size(I_total(irange))),'r'); ylabel('lux');
    ipat = find(I_total>measthresh); % indexes above measurement threshold
    
    if ~BLANKSEG % segmentation without blank frames based on timing and total number of frames
        M = params.nmeas;
        t0 = ipat(1);
        tM = ipat(end);
        interval = (tM-t0)/(M-1/3); % interval
        range = max(ceil((right-left)*interval),1); % 
        indall = round(t0+((1:M)-1/2)*interval); % mid-point [0.5]
        % meas = I(indall,:)';
        % meas = ((I(indall,:)+I(indall+1,:))/2)';
        lindall = floor(t0+((1:M)-(1-left))*interval); % left-point
        % rindall = ceil(t0+((1:M)-(1-right))*interval); % right-point
        rindall = lindall + range;
        patseg = [lindall', rindall'];
        for s = 1:M
            switch processmeth
                case {'mean','average','avr'}
                    meas(:,s) = mean(I(patseg(s,1):patseg(s,2),:),1)';
                case {'max','maximum'}
                    meas(:,s) = max(I(patseg(s,1):patseg(s,2),:),[],1)';
            end
        end
        patseg = (patseg - ipat(1))/(ipat(end) - ipat(1));
        measall = meas;
        patsegall = patseg;
    else % segmentation with blank frames based on sharp drops of blank frames
        dpat = find(diff(ipat)>2); % dramatic drop of adjacent switching
        % patno = [ipat([1;dpat+1])+50,ipat([dpat;end])-30]; % pattern indexes (from the first column to the second)
        patno = [ipat([1;dpat+1]),ipat([dpat;end])]; % pattern indexes (from the first column to the second)
        if ~includelastframe % exclude the last frame (for WinTech DMD)
            patno = patno(1:end-1,:); % exclude the last frame (excessive, 1st)
        end
        % meas(:,ich) = v(floor(0.3*patno(:,1)+0.7*patno(:,2))); % measurement for each channel
        if isfield(params, 'nmeas')
            % if known number of measurements -> total index / nmeas
            interval = (ipat(end)-ipat(1))/params.nmeas; % interval
            % patno(patno(:,2)-patno(:,1)<0.2*interval,:) = []; % exclude extra wrong data resulted from threshold
        else
            % unkown number of measurements -> estimate from average interval
            interval = round((patno(end,1)-patno(1,1)+patno(end,2)-patno(1,2))/length(patno)/2);
            patno(patno(:,2)-patno(:,1)<0.2*interval,:) = []; % exclude extra wrong data resulted from threshold
            interval = round((patno(end,1)-patno(1,1)+patno(end,2)-patno(1,2))/length(patno)/2);
        end
        if isfield(params, 'patseg')
            % turn existing patseg to relative positions
            patseg = params.patseg;
            patseg = round(patseg*(ipat(end) - ipat(1)) + ipat(1));
            nmeas  = params.nmeas;
        else
            t = 1;
            im = 1; % image segment (one if a single segment)
            patsegall = {}; % cell to include different length of each image segment
            for p = 1:length(patno)
                if p>1
                    p1 = patno(p,1);
                    p0 = patno(p-1,1); q0 = patno(p-1,2);
                    if p1-p0 > segimfactor*interval % next image segment 
                        patsegall{im} = patseg;
                        patseglen(im) = t-1;
                        im = im+1;
                        t = 1;
                    end
                    while (p1-p0>interfactor*interval) && (p1-p0<segimfactor*interval)
                        % if out of one interval range, insert a segment until within an interval
                        % disp(p1-p0);
                        p0 = p0+interval; q0 = q0+interval;
                        patseg(t,:) = [p0,q0];
                        t = t+1;
                    end
                end 
                patseg(t,:) = patno(p,:);
                t = t+1;
            end
            patsegall{im} = patseg;
            patseglen(im) = t-1;
        end

        % filter out image segment with <50% of the maximum length
        imvalid = 1;
        patseglenmax = max(patseglen);
        for im = 1:length(patsegall)
            patseg = patsegall{im};
            patsegall{im} = [];
            if patseglen(im) > 0.5*patseglenmax
                for s = 1:patseglen(im)
                    patseg(s,:) = [floor((1-left)*patseg(s,1)+left*patseg(s,2)), ceil((1-right)*patseg(s,1)+right*patseg(s,2))];
                    switch processmeth
                        case {'mean','average','avr'}
                            meas(:,s) = mean(I(patseg(s,1):patseg(s,2),:),1)';
                        case {'max','maximum'}
                            meas(:,s) = max(I(patseg(s,1):patseg(s,2),:),[],1)';
                    end
                end
                patseg = (patseg - ipat(1))/(ipat(end) - ipat(1));
                patsegall{imvalid} = patseg;
                measall{imvalid} = meas;
                imvalid = imvalid + 1;
                % NOTE: patseg, meas not reset -> a feature for same pattern sequence
            end
        end
        patsegall = patsegall(~cellfun(@isempty, patsegall));
        if length(measall) == 1
            patsegall = patsegall{1};
            measall = measall{1};
        end
    end
end

% end

