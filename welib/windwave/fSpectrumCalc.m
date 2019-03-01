function [S_average,f_average]=fSpectrumCalc(S,F1)

% Average neighboring frequency bins.

% INPUTS: Power spectrum and Nyquist Frequency
% OUTPUTS: Power spectrum and Nyquist Frequency after averaging neighboring
% frequency bins.


bins = 15; % tipically between 10 and 20
a = 10^(1/bins); % distance between neighbours
x_x=log10(F1(2)); 
fo=10^(floor(x_x));  %detection of the lowest bound of the range
F_Lower = F1(1);
F_Higher = a*fo; % higher limit of the first bin
j=0;
while (F_Lower<F1(length(F1)))
    fi=F1((F1>=F_Lower)&(F1<F_Higher));
    Si=S((F1>=F_Lower)&(F1<F_Higher));
    if (~isempty(Si))
        
        j=j+1;
        S_average(j)=mean(Si);
        f_average(j)=mean(fi);
    end
    F_Lower=F_Higher;
    F_Higher =a*F_Lower;
end

end

