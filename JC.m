function [ d2 ] = JC( X1,Xr )
%JC Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(Xr);
d2=zeros(m,1);
for i=1:m
    Xi=Xr(i,:);
    n11=0;
    n10=0;
    n01=0;
    n00=0;
    for j=1:n
        if X1(1,j)==1 && Xi(1,j)==1
            n11=n11+1;
        elseif X1(1,j)==1 && Xi(1,j)==0
            n10=n10+1;
        elseif X1(1,j)==0 && Xi(1,j)==1
            n01=n01+1;
        else
            n00=n00+1;
        end
    end
    d2(i,1)=n11/(n11+n10+n01);
end

end

