function [ x, Tx ] = zero_filling( meas, opt )
%ZERO_FILLING Reconvering sub-sampled image or image sequence signal with
%orthogonal bases by filling the non-sampled coefficients with all zeros.
%Here, we investigate three typical orthogonal bases, i.e., Walsh-Hadamard 
%transform (WHT), discrete cosine transform (DCT), and discrete Fourier 
%transform (DFT).
%  X = ZERO_FILLING(MEAS, OPT) fills the transform domain with MEAS 
%  according to the sampled indices OPT.SENSIND in full orthogonal matrix 
%  OPT.SENSMTX (if available, otherwise it will be generated according to
%  the sampling method OPT.SENSMETHOD and the corresponding forward sensing
%  matrix data generation process GENDATA.
% 
%  See also GENDATA.

rows = opt.rows;
cols = opt.cols; 
N = rows*cols; % size of the original signal
nch = size(meas, 2); % number of channels in the measurement

% zero-filling the transform domain
Tx = zeros([N nch]); % transform domain (multi-channel meas supported)
Tx(opt.sensind,:) = meas; % zero-filling

% retrieve the spatial domain using inverse transform
switch lower(opt.sensmethod)
    case {'walsh','hadamard','mixhadamard'} % Walsh-Hadamard transform (WHT)
        if ~isfield(opt, 'sensmtx')
            Wr = walsh(rows); % Walsh-Hadamard matrix in rows
            Wc = walsh(cols); % Walsh-Hadamard matrix in columns
            opt.sensmtx = kron(Wc,Wr); % full Walsh-Hadamard matrix 
        end
        x = opt.sensmtx*Tx/N;
        if strcmp(opt.sensmethod, 'mixhadamard')
            x = abs(x);
        end
    case {'dct','cosine'} % discrete cosine transform (DCT)
        if ~isfield(opt, 'sensmtx')
            Wr = dctmtx(rows); % DCT matrix in rows
            Wc = dctmtx(cols); % DCT matrix in columns
            opt.sensmtx = kron(Wc,Wr)*sqrt(rows*cols)/2; % full DCT matrix [-1,1]
        end
        % x = opt.sensmtx'*Tx;
        x = opt.sensmtx'*Tx/N*4; % need to figure out the scalar factor
    case {'fourier','dft'} % discrete Fourier transform (DFT)
        if ~isfield(opt, 'sensmtx')
            Wr = dftmtx(rows); % DFT matrix in rows
            Wc = dftmtx(cols); % DFT matrix in columns
            opt.sensmtx = kron(Wc,Wr); % full DFT matrix 
            shiftind = fftshift(reshape(1:rows*cols,[rows,cols])); % shift index of fftshift
            opt.sensmtx = opt.sensmtx(shiftind,:); % fftshift for the rows (zero frequency in the center)
        end
        x = opt.sensmtx'*Tx/N;
    case {'haar'} % Haar wavelet transform sensing matrix [-1,1] [orthogonal] 
        if ~isfield(opt, 'sensmtx')
            Wr = haarmtx(rows); % Haar wavelet transform matrix in rows
            Wc = haarmtx(cols); % Haar wavelet transform in columns
            opt.sensmtx = kron(Wc,Wr)/2; % full Haar wavelet transform matrix 
        end
        x = opt.sensmtx'*Tx/4;
    otherwise 
        error('Unsupported basis type %s!', lower(opt.sensmethod));
end

end

