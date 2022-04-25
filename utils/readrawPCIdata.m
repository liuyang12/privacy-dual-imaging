function [ volt ] = readrawPCIdata( rawPCI_data_path, params )
%READRAWPCIDATA Read raw PCI data from the capture card.
% INPUTS:
%   rawPCI_data_path - raw PCI data path
%   params           - parameters (Device | PCI8552_raw or PCI8514_raw)
% OUPUT:
%   volt             - captured voltage data, each column as a channel

fid = fopen(rawPCI_data_path);

switch params.Device
    case 'PCI8552_raw' % 12-bit capture card PCI8552
        ADBuffer = fread(fid, 'uint16');
        ADData = bitand(ADBuffer, 4095); %0x0fff=4095, lower 12 bits
        data = ((10000.0/4096) .* ADData - 5000.0);
        volt(:,1) = data(1:2:end-10); % ch1
        volt(:,2) = data(2:2:end-10); % ch2
    case 'PCI8514_raw' % 14-bit capture card PCI8514
        ADBuffer = fread(fid, 'uint16');
        ADData = bitand(ADBuffer, 16383); %0x3fff=16383, lower 14 bits
        data = ((10000.0/16384) .* ADData - 5000.0);
        volt(:,1) = data(1:4:end-20); % ch1
        volt(:,2) = data(2:4:end-20); % ch2
        volt(:,3) = data(3:4:end-20); % ch3
        volt(:,4) = data(4:4:end-20); % ch4
    otherwise
        fprintf(2,'Error! Not supported Device type!\n');   
        return;
end

fclose(fid);    % close the opened file

end

