function W = walsh(n,classname)
%WALSH Walsh matrix or Hadamard matrix in sequency order.
%   W = WALSH(N) returns is a Walsh-Hadamard matrix of order N in sequency
%   order, that is, a matrix W with elements +1 or -1 such that W'*W =
%   N*EYE(N). An N-by-N Walsh matrix with N > 2 exists only if REM(N,4)=0.
%   This function handles only the cases where N us a power of 2.
%   
%   W = WALSH(N,CLASSNAME) returns a matrix of class CLASSNAME, which can
%   be either 'single' or 'double' (the default).
%   
%   Example:
%   
%   WALSH(4) is 
%               1     1     1     1
%               1     1    -1    -1
%               1    -1    -1     1
%               1    -1     1    -1
%   See also HADAMARD.
% 
%   Copyright(C) 2017 <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>. with reference to
%   MATLAB(R) document of  <a href= 
%   "matlab:web('https://www.mathworks.com/help/signal/examples/
%   discrete-walsh-hadamard-transform.html')">Discrete Walsh-Hadamard 
%   Transform</a>
%   Last modified Dec 16, 2017.

if nargin < 2 || isempty(classname)
    classname = 'double';
end

H = hadamard(n,classname); % Hadamard matrix (in natural order)

hadIdx = 0:n-1;            % Hadamard index
m = log2(n)+1;             % Number of bits to represent the index
% Each column of the sequency index (in binary format) is given by the
% modulo-2 addition of columns of the bit-reversed Hadamard index (in
% binary format).
% Bit reversing of the binary index
binhadIdx = fliplr(dec2bin(hadIdx,m))-'0'; 
binseqIdx = zeros(n,m-1);  % Pre-allocate memory
for k = m:-1:2
    % Binary sequency index
    binseqIdx(:,k) = xor(binhadIdx(:,k),binhadIdx(:,k-1));
end
seqIdx = binseqIdx*pow2((m-1:-1:0)'); % Binary to integer sequency index
W = H(seqIdx+1,:);                    % 1-based indexing

end

