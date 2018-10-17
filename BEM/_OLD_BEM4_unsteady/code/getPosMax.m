function I =getPosMax(x)
    I=1:length(x);
    I=I(diff([x(1) x])>0 & diff([x x(end)])<0);
end