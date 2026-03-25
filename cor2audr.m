function tsf = cor2audr(z, rv, STF,BP)
% COR2AUDR complex inverse wavelet transform.
%	%	tsf = cor2audr(z, rv);
%	tsf = cor2audr(z, rv, STF, BP);
%	z: cortical representation, T-M*len(sv)-2*R, T = # of time frames
%	rv: rate vector, R-by-1
%	STF: (optional) sample  temporal frequency e.g., 125 Hz for 8 ms
%	COR2AUDR reconstructs auditory spectrum at different scales from a complex
%	cortical response Z with respect to rate vector rv.
%   The default frame resolution is 8ms i.e 125 frames/sec.
%	See also: AUD2CORS, COR2AUD , COR2AUD2

% Author: Lakshmi Krishnan (lakshmik@umd.edu), NSL, UMD
% v1.00: 19-Apr-2012

%check syntax


[T,MS,R] = size(z);
tsf = zeros(T,MS);

if (2*length(rv)) ~= R, error('size(z,3) ~= 2 * len(rv)');end;
% normalization

if nargin < 3, STF = 125; end;
K1 = length(rv);
N = T;
N1 = 2^nextpow2(N);	N2 = N1*2;

HR = zeros(2*N1, 2*K1);
for k = 1:K1
    Hr = gen_cort(rv(k), N1, STF, [k+BP K1+BP*2]);
    Hr = conj([Hr; zeros(N1, 1)]);	% SSB -> DSB
    HR(:, k+K1) = Hr;               % downward
    HR(:, k) = [Hr(1); conj(flipud(Hr(2:N2)))];	% upward
    HR(N1+1, k) = abs(HR(N1+2, k));
end

Z_cum = zeros(N2,MS);
Hsum = zeros(N2,1);
% cumulation



for r = 1:length(rv),
    
    for sgn = [-1 1],
        
        % rate filtering modification
        if sgn < 0,
            
            for i = 1:size(z,2)
                ZF	= fft(squeeze(z(:,i,r+length(rv))),N2);
                Z_cum(:,i) = Z_cum(:,i) + ZF.*HR(:,r+K1);
                
            end
            Hsum       = Hsum +HR(:,r+K1).^2;
        else
            
            
            for i = 1:size(z,2)
                ZF	= fft(squeeze(z(:,i,r)), N2);
                Z_cum(:,i) = Z_cum(:,i) + ZF.*HR(:,r);
                
            end
           Hsum       = Hsum +HR(:,r).^2;
        end
        
        
    end
end


NORM =0.99;
Hsum(1,:) = Hsum(1,:)*2;
sumH = sum(Hsum,1);
Hsum = NORM * Hsum + (1 - NORM) * max(Hsum(:));
Hsum = Hsum ./ sum(Hsum(:)) * sumH + 1;	 % normalization for DC

for i = 1:size(z,2)
    Z_cum(:,i) = Z_cum(:,i) ./ Hsum;
    tmp =  ifft(Z_cum(:,i));
    tsf(:,i) = tmp(1:T);
    
end


