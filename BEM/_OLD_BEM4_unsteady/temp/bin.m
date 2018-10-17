function [X,Y]=bin(x,y,n)
    X=zeros(1,n);
    Y=zeros(1,n);
    step=(max(x)-min(x))/n;
    for i=1:n
        X(i)=mean(x( x>=(i-1)*step+min(x) & x<i*step+min(x)  ));
        Y(i)=mean(y(  x>=(i-1)*step+min(x) & x<i*step+min(x)   ));
    end
end